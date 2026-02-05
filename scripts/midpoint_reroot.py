"""
Use BTE to read in an UShER protobuf file and traverse it to find the midpoint node.
Run the matUtils command to reroot to that node, producing a modified reference fasta for the new root sequence.
Run usher_to_taxonium to visualize the rerooted tree.
"""

import argparse
import bte
import gzip
import os
import subprocess
import sys
import tempfile
import alter_gbff


def farthest_descendant_excluding(node: bte.MATNode, current_distance: int, current_path: list, max_distance: int, longest_path: list,
                                  include_leaves: bool, excluded_node: bte.MATNode):
    """Recursively descend from node to find the path to the descendant with the greatest sum of mutation counts relative to node,
    but skipping excluded_node."""
    for child in node.children:
        if child.id != excluded_node.id:
            if len(child.children) > 0 or include_leaves:
                child_distance = current_distance + len(child.mutations)
                child_path = current_path + [child]
                if child_distance > max_distance:
                    max_distance = child_distance
                    longest_path = child_path
                max_distance, longest_path = farthest_descendant_excluding(child, child_distance, child_path, max_distance, longest_path,
                                                                           include_leaves, excluded_node)
    return max_distance, longest_path


def farthest_descendant(node: bte.MATNode, current_distance: int, current_path: list, max_distance: int, longest_path: list, include_leaves: bool):
    """Recursively descend from node to find the path to the descendant with the greatest sum of mutation counts relative to node."""
    # For the excluded node, just pass in node itself because that will never match a descendant.
    return farthest_descendant_excluding(node, current_distance, current_path, max_distance, longest_path, include_leaves, node)


def find_farthest_from_node(source_node: bte.MATNode, include_leaves: bool):
    """Find the farthest descendant of source_node, then search outward through the rest of the tree in case a non-descendant is farther."""
    max_distance, longest_path = farthest_descendant(source_node, 0, [source_node], 0, [source_node], include_leaves)
    print(f"Farthest descendant from {source_node.id} is {longest_path[-1].id} with a distance of {max_distance}")
    search_node = source_node
    current_distance = 0
    current_path = [source_node]
    while search_node.parent is not None:
        current_distance += len(search_node.mutations)
        exclude_node = search_node
        search_node = search_node.parent
        current_path.append(search_node)
        max_distance, longest_path = farthest_descendant_excluding(search_node, current_distance, current_path, max_distance, longest_path,
                                                                   include_leaves, exclude_node)
    print(f"Farthest node from {source_node.id} is {longest_path[-1].id} with a distance of {max_distance}")
    return max_distance, longest_path


def find_path_midpoint(path: list):
    """Return the midpoint of path, i.e. the parent node whose child edge distance spans the midpoint of the total path length."""
    distances = []
    path_length = 0
    for ix in range(len(path) - 1):
        if path[ix+1].parent and path[ix].id == path[ix+1].parent.id:
            distance = len(path[ix+1].mutations)
        elif path[ix].parent and path[ix+1].id == path[ix].parent.id:
            distance = len(path[ix].mutations)
        else:
            print(f"Consecutive nodes [{ix}]={path[ix].id} and [{ix+1}]={path[ix+1].id} do not have a parent-child relationship. {path[ix].id}'s parent is {path[ix].parent.id} and {path[ix+1].id}'s parent is {path[ix+1].parent.id}", file=sys.stderr)
            sys.exit(1)
        distances.append(distance)
        path_length += distance
    midpoint = path_length / 2.0
    print(f"midpoint distance is {midpoint}")
    midpoint_ix = 0
    current_distance = 0
    for ix in range(len(distances)):
        current_distance += distances[ix]
        if current_distance >= midpoint:
            midpoint_ix = ix
            break
    if path[midpoint_ix] == path[midpoint_ix+1].parent:
        print(f"midpoint is {path[midpoint_ix].id} (child {path[midpoint_ix+1].id}) with distance {current_distance}")
        return path[midpoint_ix]
    else:
        print(f"midpoint is {path[midpoint_ix+1].id} (child {path[midpoint_ix].id}) with distance {current_distance}")
        return path[midpoint_ix+1]


def find_tree_midpoint(input_tree_file: str, include_leaves: bool):
    """Find the midpoint node of an UShER tree by reading it in with BTE, finding the node that is farthest from root,
    finding the node that is farthest from that node, and taking the node of the midpoint of that path.  Return the node's ID."""
    tree = bte.MATree(input_tree_file)
    max_distance_root, longest_path_root = find_farthest_from_node(tree.root, include_leaves)
    farthest_from_root = longest_path_root[-1]
    max_distance, longest_path = find_farthest_from_node(farthest_from_root, include_leaves)
    if max_distance < max_distance_root:
        print(f"Note: max_distance {max_distance} skipping root-child is less than max distance from root {max_distance_root}; using max from root")
        longest_path = longest_path_root
    return find_path_midpoint(longest_path).id


def run_command(command: list):
    """Run a command (list of arg string), exit on failure."""
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"command ({" ".join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)


def reroot_tree(input_tree_file: str, input_ref_fasta: str, new_root_node: str, output_base: str):
    """Run matUtils to reroot the tree and produce modified reference fasta."""
    output_tree_file = output_base + ".pb.gz"
    output_fasta_file = output_base + ".fasta"
    command = ["/cluster/home/angie/github/usher/build/matUtils", "extract",
               "--input-mat", input_tree_file,
               "--input-fasta", input_ref_fasta,
               "--reroot", new_root_node,
               "--write-reroot-reference", output_fasta_file,
               "--write-mat", output_tree_file
               ]
    run_command(command)
    return output_tree_file, output_fasta_file


def get_accession(fasta_file: str):
    """Expect the first line to begin with > and return the first space-separated word after ">"."""
    with open(fasta_file) as f:
        line = f.readline().strip()
        if line.startswith(">"):
            accession = line[1:].split(" ")[0]
            return accession
        else:
            print(f"Unable to parse ID from first line of fasta file {fasta_file}:\n{line}", file=sys.stderr)
            sys.exit(1)


def make_gbff(gbff_file: str, accession: str, fasta_file: str, output_base: str):
    """Modify the original reference GBFF file using the new rerooted reference sequence."""
    output_gbff_file = output_base + ".gbff"
    alter_gbff.alter_gbff_file(gbff_file, accession, fasta_file, output_gbff_file)
    return output_gbff_file


def tweak_metadata(metadata_file: str):
    """Strip the long-names first column from the metadata file so the accession column is first."""
    # Make a temporary file that will need to be deleted later
    tweaked_metadata_file = ""
    with tempfile.NamedTemporaryFile(suffix='.tsv.gz', delete=False) as tmp:
        tweaked_metadata_file = tmp.name
        with gzip.GzipFile(fileobj=tmp, mode='w') as gz_file:
            with gzip.open(metadata_file, 'rt') if metadata_file.endswith('.gz') else open(metadata_file) as f:
                for line in f:
                    words = line.split("\t")
                    gz_file.write("\t".join(words[1:]).encode('utf-8'))
    return tweaked_metadata_file


def get_header_words(tsv_file):
    """Return a list of the words in a TSV file header"""
    with gzip.open(tsv_file, 'rt') if tsv_file.endswith('.gz') else open(tsv_file) as f:
        line = f.readline().rstrip()
        return line.split("\t")


def make_taxonium(tree_file: str, metadata_file: str, gbff_file: str, output_base: str):
    """Run usher_to_taxonium so we can view the rerooted tree"""
    output_taxonium_file = output_base + ".jsonl.gz"
    modified_metadata_file = tweak_metadata(metadata_file)
    columns = get_header_words(modified_metadata_file)
    command = ["usher_to_taxonium",
               "--input", tree_file,
               "--metadata", modified_metadata_file,
               "--key_column", columns[0],
               "--columns", ",".join(columns),
               "--genbank", gbff_file,
               "--name_internal_nodes",
               "--title", output_base,
               "--output", output_taxonium_file,
               ]
    run_command(command)
    os.remove(modified_metadata_file)


def main():
    parser = argparse.ArgumentParser(description="Find midpoint of an UShER tree.  Make a rerooted tree, fasta, GBFF, and taxonium file.")
    parser.add_argument("-i", "--input_tree", required=True,
                        help="Input UShER protobuf tree file")
    parser.add_argument("-g", "--ref_gbff", required=True,
                        help="Input GenBank Flat File describing reference for the tree")
    parser.add_argument("-m", "--metadata", required=True,
                        help="TSV metadata file for tree.  First column must match tree leaf names.")
    parser.add_argument("-r", "--ref_fasta", required=True,
                        help="Input reference fasta file for tree")
    parser.add_argument("-o", "--output_base", required=True,
                        help="Output file base name for .pb.gz, .fasta and .jsonl.gz output files.")
    parser.add_argument("--include_leaves", action='store_true',
                        help="Include leaf distances (private mutations), instead of considering only internal nodes (the default)")
    args = parser.parse_args()

    midpoint_node = find_tree_midpoint(args.input_tree, args.include_leaves)
    rerooted_tree, modified_ref_fasta = reroot_tree(args.input_tree, args.ref_fasta, midpoint_node, args.output_base)
    ref_accession = get_accession(args.ref_fasta)
    modified_gbff = make_gbff(args.ref_gbff, ref_accession, modified_ref_fasta, args.output_base)
    make_taxonium(rerooted_tree, args.metadata, modified_gbff, args.output_base)


if __name__ == "__main__":
    main()
