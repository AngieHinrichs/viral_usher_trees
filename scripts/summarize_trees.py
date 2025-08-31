#!/usr/bin/env python3
"""
Given a metadata TSV file describing RefSeqs with tree names, find subdirectories of trees/ that have tree files
and output_stats.tsv files, and write a TSV file that adds the size of each tree (in number of tips) to metadata.
"""

import argparse
import csv
import os
import sys

import viral_usher_trees


metadata_required_columns = ["accession", "tree_name", "organism", "isolate", "strain", "serotype", "segment", "Taxonomy_ID"]


def get_tip_count(output_stats_tsv):
    """Return value in tree_tip_count column (should be only one data row); return -1 if file is not as expected"""
    tip_count = -1
    with open(output_stats_tsv, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        if "tree_tip_count" not in fieldnames:
            return -1
        for row in reader:
            tip_count = int(row["tree_tip_count"])
    return tip_count


def summarize_trees(input_tsv, output_tsv):
    """Use tree metadata TSV and trees/*/output_stats.tsv to make summary TSV"""
    with open(input_tsv, 'r', newline='', encoding='utf-8') as f, open(output_tsv, "w") as outf:
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
        writer = csv.writer(outf, delimiter='\t', lineterminator='\n')
        out_header = metadata_required_columns.copy()
        out_header.append("tip_count")
        writer.writerow(out_header)
        for row in reader:
            acc = row["accession"]
            tree = row["tree_name"]
            org = row["organism"]
            iso = row["isolate"]
            stn = row["strain"]
            ser = row["serotype"]
            seg = row["segment"]
            tax = row["Taxonomy_ID"]
            subdir_path = viral_usher_trees.trees_dir + "/" + tree
            if not os.path.isdir(subdir_path):
                continue
            output_stats_path = subdir_path + "/" + "output_stats.tsv"
            if not os.path.exists(output_stats_path):
                continue
            tip_count = get_tip_count(output_stats_path)
            if tip_count < 0:
                print(f"Contents of {output_stats_path} are not as expected!  Skipping {tree}.", file=sys.stderr)
                continue

            writer.writerow([acc, tree, org, iso, stn, ser, seg, tax, tip_count])


def main():
    parser = argparse.ArgumentParser(description="Use tree metadata TSV and trees/*/output_stats.tsv to make summary TSV")
    parser.add_argument('-i', '--input', required=True,
                        help="Input RefSeq tree metadata TSV file, must include columns " + ", ".join(metadata_required_columns))
    parser.add_argument('-o', '--output', required=True,
                        help="Output tree summary TSV file, columns are required input columns and tip_count")
    args = parser.parse_args()
    viral_usher_trees.check_top_level_dir()
    summarize_trees(args.input, args.output)


if __name__ == "__main__":
    main()
