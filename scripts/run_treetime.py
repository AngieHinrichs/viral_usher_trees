#!/usr/bin/env python3
"""
If the proportion of non-empty date values in metadata.tsv.gz is high enough,
then get the RefSeq length from output_stats.tsv, make a newick file with
branch lengths scaled as substitutions per site, make a dates.csv file with
added -XX for missing month and/or day, and run treetime to reroot the tree.
Then apply the same rerooting to optimized.pb.gz.
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
import treeswift

import viral_usher_trees
import alter_gbff

default_min_real_dates = 0.8


def get_dates(subdir_path):
    """Scan metadata.tsv.gz, return a dict of names -> dates and the proportion of real date values to total count"""
    name_to_date = {}
    with gzip.open(subdir_path + "/metadata.tsv.gz", "rt") as f:
        header = [field.strip() for field in f.readline().split('\t')]
        name_idx = header.index('accession')
        if name_idx < 0:
            print(subdir_path + "/metadata.tsv.gz does not have accession column", file=sys.stderr)
            sys.exit(1)
        date_idx = header.index('date')
        if date_idx < 0:
            print(subdir_path + "/metadata.tsv.gz does not have date column", file=sys.stderr)
            sys.exit(1)
        line_count = 0
        real_date_count = 0
        for line in f:
            fields = line.rstrip().split("\t")
            name = fields[name_idx]
            date = fields[date_idx]
            if re.match('^[0-9]{4}', date):
                name_to_date[name] = date
                real_date_count += 1
            line_count += 1
    return name_to_date, (real_date_count / line_count)


def scale_branch_lengths(subdir_path, newick_out, refseq_len):
    """Use treeswift (included with taxoniumtools) to scale branch lengths as expected by treetime"""
    tmp_newick_path = subdir_path + "/optimized.nwk"
    command = ["matUtils", "extract",
               "-i", subdir_path + "/optimized.pb.gz",
               "-t", tmp_newick_path]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"matUtils command ({" ".join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)
    tree = treeswift.read_tree_newick(tmp_newick_path)
    tree.scale_edges(1.0 / refseq_len)
    tree.write_tree_newick(subdir_path + "/" + newick_out)
    os.remove(tmp_newick_path)


def make_dates_csv(subdir_path, dates_out, name_to_date):
    """Write a CSV file of names and dates, fillin in missing month and/or day with -XX"""
    with open(subdir_path + "/" + dates_out, "w") as f:
        f.write(",".join(["name", "date"]) + "\n")
        for name, date in name_to_date.items():
            if re.match('^[0-9]{4}$', date):
                date += "-XX-XX"
            elif re.match('^[0-9]{4}-[0-9]{2}$', date):
                date += "-XX"
            f.write(",".join([name, date]) + "\n")


def get_refseq_len(subdir_path):
    """Get refseq_length value from output_stats.tsv"""
    with open(subdir_path + "/output_stats.tsv", "r", newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        if "ref_length" not in fieldnames:
            print("output_stats.tsv file does not have ref_length column", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            ref_length = int(row["ref_length"])
    return ref_length


def run_treetime(tree, min_real_dates):
    """Format input for treetime clock and run it"""
    subdir_path = viral_usher_trees.trees_dir + "/" + tree
    if not os.path.isdir(subdir_path):
        print(f"{viral_usher_trees.trees_dir} does not have a subdirectory {tree}, check spelling", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(subdir_path + "/metadata.tsv.gz") or not os.path.exists(subdir_path + "/output_stats.tsv"):
        print(f"Expected files not found in {subdir_path}, has the tree been built?", file=sys.stderr)
        sys.exit(1)
    name_to_date, real_dates_proportion = get_dates(subdir_path)
    if real_dates_proportion < min_real_dates:
        print(f"Tree {tree} has too low a proportion of dates ({real_dates_proportion:.2f} < {min_real_dates:.2f}), not running treetime. ")
        sys.exit(0)
    refseq_len = get_refseq_len(subdir_path)
    scale_branch_lengths(subdir_path, "optimized.scaled.nwk", refseq_len)
    make_dates_csv(subdir_path, "dates.csv", name_to_date)
    command = ["treetime", "clock",
               "--sequence-length", str(refseq_len),
               "--tree",  subdir_path + "/optimized.scaled.nwk",
               "--dates", subdir_path + "/dates.csv",
               "--outdir", subdir_path + "/treetime_out"]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"treetime command ({" ".join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)


def get_oldest_node(tree):
    """Extract the oldest internal node date from treetime's output file rtt.csv"""
    # oldest_node=$(grep ^node_ treetime_out/rtt.csv | sort -t, -k2n | head -1 | cut -d, -f 1)
    rtt_csv_path = "/".join([viral_usher_trees.trees_dir, tree, "treetime_out", "rtt.csv"])
    min_date = None
    oldest_node = None
    with open(rtt_csv_path) as f:
        for line in f:
            columns = line.split(",")
            if columns[0].startswith("node_"):
                date = float(columns[1].strip())
                if min_date is None or date < min_date:
                    min_date = date
                    oldest_node = columns[0].strip()
    if oldest_node is None:
        print(f"Failed to get oldest node from {rtt_csv_path}", file=sys.stderr)
        sys.exit(1)
    return oldest_node


def get_refseq_acc(config_path):
    """Return the value of refseq_acc from config.toml"""
    refseq_acc = None
    with open(config_path) as f:
        for line in f:
            if line.startswith("refseq_acc = "):
                refseq_acc = line.split("'")[1]
    if refseq_acc is None:
        print(f"Failed to find refseq_acc in {config_path}", file=sys.stderr)
        sys.exit(1)
    return refseq_acc


def reroot_tree(tree, oldest_node):
    """Reroot optimized.pb.gz to oldest_node using matUtils and return path to rerooted .pb.gz."""
    input_path = "/".join([viral_usher_trees.trees_dir, tree, "optimized.pb.gz"])
    rerooted_tree_path = "/".join([viral_usher_trees.trees_dir, tree, "timetree_rerooted.pb.gz"])
    config_path = "/".join([viral_usher_trees.trees_dir, tree, "config.toml"])
    refseq_acc = get_refseq_acc(config_path)
    ref_path = "/".join([viral_usher_trees.trees_dir, tree, refseq_acc + ".fa"])
    modified_ref_path = "/".join([viral_usher_trees.trees_dir, tree, "treetime_rerooted_" + refseq_acc + ".fa"])
    command = ["matUtils", "extract",
               "-i", input_path,
               "--reroot", oldest_node,
               "--input-fasta", ref_path,
               "--write-reroot-reference", modified_ref_path,
               "-o", rerooted_tree_path]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"matUtils command ({" ".join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)
    gbff_path = "/".join([viral_usher_trees.trees_dir, tree, refseq_acc + ".gbff"])
    rerooted_gbff_path = "/".join([viral_usher_trees.trees_dir, tree, "treetime_rerooted_" + refseq_acc + ".gbff"])
    alter_gbff.alter_gbff_file(gbff_path, refseq_acc, ref_path, rerooted_gbff_path)
    return rerooted_tree_path, rerooted_gbff_path


def get_columns_string(tsv_path):
    """Return a comma-separated string listing columns of TSV file."""
    # columns=$(gunzip -c metadata.tsv.gz | head -1 | sed -re 's/\t/,/g')
    with gzip.open(tsv_path, "rt") if tsv_path.endswith(".gz") else open(tsv_path) as f:
        header = f.readline().split('\t')
        columns = ",".join([col.strip() for col in header])
        return columns
    print(f"Failed to get columns from {tsv_path}", file=sys.stderr)
    sys.exit(1)


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


def make_taxonium(tree, rerooted_tree_path, rerooted_gbff_path):
    """Run usher_to_taxonium on rerooted_tree"""
    input_metadata_path = "/".join([viral_usher_trees.trees_dir, tree, "metadata.tsv.gz"])
    tmp_metadata_path = tweak_metadata(input_metadata_path)
    columns = get_columns_string(tmp_metadata_path)
    taxonium_jsonl_path = "/".join([viral_usher_trees.trees_dir, tree, "timetree_rerooted.jsonl.gz"])
    command = ["usher_to_taxonium",
               "--input", rerooted_tree_path,
               "--metadata", tmp_metadata_path,
               "--key_column", "accession",
               "--genbank", rerooted_gbff_path,
               "--columns", columns,
               "--name_internal_nodes",
               "--title", "Treetime-rerooted " + tree,
               "--output", taxonium_jsonl_path]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"usher_to_taxonium command ({" ".join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)
    os.remove(tmp_metadata_path)


def main():
    parser = argparse.ArgumentParser(description="Use tree metadata TSV and trees/*/output_stats.tsv to make summary TSV")
    parser.add_argument('-t', '--tree', required=True,
                        help="Tree name (must match a name from the tree_name column of tree_metadata.tsv)")
    parser.add_argument('-m', '--min_real_dates',
                        help=f"Minimum proportion of dates in metadata.tsv.gz that have real values (default: {default_min_real_dates})")
    args = parser.parse_args()
    viral_usher_trees.check_top_level_dir()
    min_real_dates = args.min_real_dates if args.min_real_dates else default_min_real_dates
    run_treetime(args.tree, min_real_dates)
    oldest_node = get_oldest_node(args.tree)
    rerooted_tree, rerooted_gbff = reroot_tree(args.tree, oldest_node)
    make_taxonium(args.tree, rerooted_tree, rerooted_gbff)


if __name__ == "__main__":
    main()
