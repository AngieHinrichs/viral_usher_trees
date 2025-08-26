#!/usr/bin/env python3
"""
Make tree subdirectories and config files (if they don't already exist) for viruses in input TSV.
"""

import argparse
import csv
import os
import subprocess
import sys

metadata_required_columns = ["accession", "tree_name", "assembly", "Taxonomy_ID"]
config_name = "config.toml"


def generate_config(subdir_path: str, tree_name: str, refseq_acc: str, refseq_assembly: str, taxid: str) -> bool:
    """Make subdir_path/config.toml using viral_usher init with command-line args
    to skip the interactive process."""
    command = ["viral_usher", "init",
               "--refseq", refseq_acc,
               "--taxonomy_id", taxid,
               # Nextclade dataset search is not working well enough, and command fails if there are multiple matches,
               # so just say no to nextclade by default; add back manually to config.toml files where applicable.
               "--nextclade_dataset", "",
               "--workdir", subdir_path,
               "--config", subdir_path + "/" + config_name]
    print(" ".join(command))
    try:
        subprocess.run(command, check=True)
        return True
    except FileNotFoundError as e:
        print(f"viral_usher needs to be in your path. ({e})", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError:
        return False


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
            subdir_path = "trees/" + tree
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
            config_path = subdir_path + "/" + config_name
            if not os.path.exists(config_path):
                if not generate_config(subdir_path, tree, acc, asm, taxid):
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
    if not os.path.isdir("trees"):
        print("Can't find trees directory.  This script must be run in top level of repo.", file=sys.stderr)
        sys.exit(1)
    generate_configs(args.input)


if __name__ == "__main__":
    main()
