#!/usr/bin/env python3
"""
Make tree subdirectories and config files (if they don't already exist) for viruses in input TSV.
"""

import argparse
import csv
import os
import sys

import viral_usher_trees

metadata_required_columns = ["accession", "tree_name", "assembly", "Taxonomy_ID"]


def generate_configs(input_tsv: str):
    """For each tree specification in input_tsv, if there isn't already a subdir with a config file, make one."""
    with open(input_tsv, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        missing = False
        for required_col in metadata_required_columns:
            if required_col not in fieldnames:
                missing = True
                print(f"Input file {input_tsv} must have a column named '{required_col}'.", file=sys.stderr)
        if missing:
            print("Sorry, can't proceed unless all required columns are present", file=sys.stderr)
            sys.exit(1)
        failures = []
        for row in reader:
            acc = row["accession"]
            tree = row["tree_name"]
            asm = row["assembly"]
            taxid = row["Taxonomy_ID"]
            subdir_path = viral_usher_trees.trees_dir + "/" + tree
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
            config_path = subdir_path + "/" + viral_usher_trees.config_name
            if not os.path.exists(config_path):
                if not viral_usher_trees.generate_config(subdir_path, tree, acc, asm, taxid):
                    failures.append(tree)
        if failures:
            print("\n*** Some viral_usher init commands failed:", file=sys.stderr)
            for tree in failures:
                print("  " + tree, file=sys.stderr)
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Make tree subdirectories and config files driven by tree metadata TSV")
    parser.add_argument('-i', '--input', required=True,
                        help="Input RefSeq tree metadata TSV file, must include columns " + ", ".join(metadata_required_columns))
    args = parser.parse_args()
    viral_usher_trees.check_top_level_dir()
    generate_configs(args.input)


if __name__ == "__main__":
    main()
