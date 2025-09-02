#!/usr/bin/env python3
"""
Read tree_metadata.tsv to find trees with fewer than a minimum number of sequences.
Add their RefSeq accessions to rank_refseqs/exclude_accessions.tsv.  git rm their files
in the appropriate trees/ subdirectory.
"""

import argparse
import csv
import subprocess
import sys

import viral_usher_trees

default_min_sequences = 100
exclude_file = "rank_refseqs/exclude_accessions.tsv"
metadata_required_columns = ["accession", "tree_name", "tip_count"]


def remove_tree_files(tree):
    """git rm all files in the trees/{tree} subdirectory"""
    subdir_path = viral_usher_trees.trees_dir + "/" + tree
    command = ["git", "rm", f"{subdir_path}/*"]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as e:
        print(f"git needs to be in your path. ({e})", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print("Command '" + " ".join(command) + f"' failed:\n{e}", file=sys.stderr)
        sys.exit(1)


def exclude_small_trees(input_tsv: str, min_sequences: int):
    """Read metadata to find trees with < min_sequences.  Add accession to exclude file, remove tree files."""
    with open(input_tsv, 'r', newline='', encoding='utf-8') as f, open(exclude_file, 'a') as exf:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        viral_usher_trees.check_required_columns(fieldnames, metadata_required_columns, input_tsv)
        exclude_count = 0
        for row in reader:
            acc = row["accession"]
            tree = row["tree_name"]
            tip_count = int(row["tip_count"])
            if tip_count < min_sequences:
                print(f"Excluding and removing {tree} with {tip_count} sequences", file=sys.stderr)
                exf.write(f"{acc}\tFinal tree has too few sequences ({tip_count})\n")
                remove_tree_files(tree)
                exclude_count += 1
        print(f"Excluded and removed {exclude_count} trees", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Exclude & remove trees with too few sequences")
    parser.add_argument('-i', '--input', required=True,
                        help="Input tree metadata TSV file, must include columns " + ", ".join(metadata_required_columns))
    parser.add_argument('-m', '--min_sequences', type=int,
                        help=f"Minimum number of sequences in order to keep a tree (default: {default_min_sequences})")
    args = parser.parse_args()
    viral_usher_trees.check_top_level_dir()
    min_sequences = args.min_sequences if args.min_sequences is not None else default_min_sequences
    exclude_small_trees(args.input, min_sequences)


if __name__ == "__main__":
    main()
