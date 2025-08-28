#!/usr/bin/env python3
"""
Given a metadata file that describes trees, assumed to be ordered roughly largest-first, and a number of concurrent jobs,
make one file per job that lists the trees to build or update.  If a tree doesn't yet have a config file then make one.
"""

import argparse
import csv
import os
import sys

import viral_usher_trees


metadata_required_columns = ["accession", "tree_name", "assembly", "Taxonomy_ID"]


def partition_by_index_mod_count(lst, count):
    """Return a list of count lists, such that [i] has elements of lst whose index mod count is i.
    If lst is sorted and its length >> count, the resulting partitions will have similar distributions,
    at least more similar than if the list were linearly chopped into count segments."""
    partitions = [[] for _ in range(count)]
    for i, item in enumerate(lst):
        partitions[i % count].append(item)
    return partitions


def partition_jobs(input_tsv, job_count):
    """Partition trees in input_tsv into job_count lists/output files.  Generate configs if missing."""
    trees_to_build = []
    config_failures = []
    with open(input_tsv, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        missing = False
        for required_col in metadata_required_columns:
            if required_col not in fieldnames:
                missing = True
                print(f"Input file {input_tsv} must have a column named '{required_col}'.", file=sys.stderr)
        if missing:
            print(f"Sorry, can't proceed unless all required columns are present in {input_tsv}", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            acc = row["accession"]
            tree = row["tree_name"]
            asm = row["assembly"]
            taxid = row["Taxonomy_ID"]
            subdir_path = viral_usher_trees.trees_dir + "/" + tree
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
            config_path = subdir_path + "/" + viral_usher_trees.config_name
            if os.path.exists(config_path) or viral_usher_trees.generate_config(subdir_path, tree, acc, asm, taxid):
                trees_to_build.append(tree)
            else:
                config_failures.append(tree)
    # Config failures require investigation; meanwhile we want to continue with the rest of the build.
    # Write a file of trees for which we couldn't generate a config:
    with open("config_failures.txt", "w") as fail_out:
        for tree in config_failures:
            fail_out.write(tree + "\n")
    # Partition the trees with config files into job_count batches and write job list files.
    partitions = partition_by_index_mod_count(trees_to_build, job_count)
    for i in range(job_count):
        with open(viral_usher_trees.job_file_name(i), "w") as job_out:
            for tree in partitions[i]:
                job_out.write(tree + "\n")


def main():
    parser = argparse.ArgumentParser(description="Partition trees for concurrent jobs, driven by tree metadata TSV.  " +
                                     "If a tree doesn't already have a config file then this attempts to generate one.  " +
                                     "This creates files config_failures.txt and job_list.*.txt in the current directory.")
    parser.add_argument('-i', '--input', required=True,
                        help="Input RefSeq tree metadata TSV file, must include columns " + ", ".join(metadata_required_columns))
    parser.add_argument('-j', '--jobs', type=int, required=True,
                        help="Number of concurrent jobs.  Output files job_list.[0..jobs-1].txt will be created.")
    args = parser.parse_args()
    viral_usher_trees.check_top_level_dir()
    partition_jobs(args.input, args.jobs)


if __name__ == "__main__":
    main()
