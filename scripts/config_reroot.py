#!/usr/bin/env python3
"""
Update the given tree's config.toml as follows:
* comment out refseq_acc so there's a record of what RefSeq was used to make the tree
* comment out refseq_assembly
* add ref_fasta = given rerooted-reference fasta file
* add ref_gbff = given rerooted-reference GenBank file
* update viral_usher_version to 0.10.0 because we're using new keywords ref_fasta and ref_gbff
"""

import argparse
import os
import re
import sys

import viral_usher_trees


def config_reroot(tree:str, ref_fasta:str, ref_gbff:str) -> str:
    lines_out = []
    config_path = "/".join([viral_usher_trees.trees_dir, tree, "config.toml"])
    with open(config_path) as cfg_in:
        for line in cfg_in:
            m = re.match(r"^([A-Za-z0-9_]+) = '[^']*'", line)
            if m:
                setting = m.groups()[0]
                if setting == 'viral_usher_version':
                    line_out = f"{setting} = '0.10.0'\n"
                elif setting == 'refseq_acc' or setting == 'refseq_assembly':
                    line_out = "#" + line
                    if setting == 'refseq_assembly':
                        line_out += f"ref_fasta = '{ref_fasta}'\n"
                        line_out += f"ref_gbff = '{ref_gbff}'\n"
                else:
                    line_out = line
                lines_out.append(line_out)
            else:
                lines_out.append(line)
    contents_out = "".join(lines_out)
    with open(config_path, "w") as cfg_out:
        print(contents_out, file=cfg_out, end="")


def main():
    parser = argparse.ArgumentParser(description="Update config.toml to use ref_fasta and ref_gbff instead of refseq_acc and refseq_assembly")
    parser.add_argument('-t', '--tree', required=True,
                        help="Tree name (must match a name from the tree_name column of tree_metadata.tsv)")
    parser.add_argument('-f', '--ref_fasta', required=True,
                        help="Rerooted-reference fasta file")
    parser.add_argument('-g', '--ref_gbff', required=True,
                        help="Rerooted-reference GenBank file")
    args = parser.parse_args()
    # Make sure we are in the top-level directory
    viral_usher_trees.check_top_level_dir()
    # Make sure ref_fasta and ref_gbff are paths that exist in the top-level directory
    # and start with "./" because build_one.sh expects that.
    for f in [args.ref_fasta, args.ref_gbff]:
        if not os.path.exists(f):
            print(f"File '{f}' not found in top-level directory", file=sys.stderr)
            sys.exit(1)
        if not f.startswith("./"):
            print(f"File '{f}' path needs to start with './'", file=sys.stderr)
            sys.exit(1)
    config_reroot(args.tree, args.ref_fasta, args.ref_gbff)


if __name__ == "__main__":
    main()
