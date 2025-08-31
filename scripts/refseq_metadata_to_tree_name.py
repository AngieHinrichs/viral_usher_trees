#!/usr/bin/env python3
"""
Make a meaningful-as-possible name for each RefSeq based on available metadata,
output TSV mapping of RefSeq accession to name.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from typing import NamedTuple

metadata_required_columns = ["accession", "organism", "assembly", "isolate", "strain", "serotype", "segment"]


class Sequence(NamedTuple):
    accession: str
    segment: str


class Assembly(NamedTuple):
    sequences: list
    isolates: list
    strains: list
    serotypes: list


def read_metadata_tsv(input_file: str) -> defaultdict[str, defaultdict[str, Assembly]]:
    """Return a dictionary that maps each organism to a dictionary that maps assembly identifiers to
    Assembly structures containing one or more Sequence tuples and various distinguishing attributes."""
    organism_assemblies = defaultdict(lambda: defaultdict(lambda: Assembly(sequences=[], isolates=[], strains=[], serotypes=[])))
    with open(input_file, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        missing = False
        for required_col in metadata_required_columns:
            if required_col not in fieldnames:
                missing = True
                print(f"Metadata file {input_file} must have a column named '{required_col}'.", file=sys.stderr)
        if missing:
            print("Sorry, can't proceed unless all required columns are present", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            acc = row["accession"]
            org = row["organism"]
            asm = row["assembly"]
            iso = row["isolate"]
            stn = row["strain"]
            ser = row["serotype"]
            seg = row["segment"]
            organism_assemblies[org][asm].sequences.append(Sequence(accession=acc, segment=seg))
            if iso and iso not in organism_assemblies[org][asm].isolates:
                organism_assemblies[org][asm].isolates.append(iso)
            if stn and stn not in organism_assemblies[org][asm].strains:
                organism_assemblies[org][asm].strains.append(stn)
            if ser and ser not in organism_assemblies[org][asm].serotypes:
                organism_assemblies[org][asm].serotypes.append(ser)
    return organism_assemblies


def sanitize_name(name):
    """Replace characters with special meaning in Newick or filesystems/directories with _"""
    for char in " ()[]{}-.:;/'"'"':
        name = name.replace(char, '_')
    name = re.sub(r'__+', '_', name)
    name = re.sub(r'^_', '', name)
    name = re.sub(r'_$', '', name)
    return name


def choose_assembly_suffix(assembly):
    """Pick the shortest available distinguishing attribute, sanitize & return it"""
    if assembly.serotypes:
        suffix = min(assembly.serotypes, key=len)
    elif assembly.isolates:
        suffix = min(assembly.isolates, key=len)
    elif assembly.strains:
        suffix = min(assembly.strains, key=len)
    else:
        suffix = ""
    return sanitize_name(suffix)


def main():
    parser = argparse.ArgumentParser(
        description="Make unique names for RefSeqs derived from metadata")
    parser.add_argument('-m', '--metadata', required=True,
                        help="Input RefSeq metadata TSV file, must have columns " + ", ".join(metadata_required_columns))
    parser.add_argument('-o', '--output', required=True,
                        help="Output TSV mapping RefSeq accession to name")
    args = parser.parse_args()

    # First collate RefSeqs by organism and assembly.
    organism_assemblies = read_metadata_tsv(args.metadata)

    # Now for each organism, if there are multiple assemblies and/or multiple segments, find available suffixes to add
    # until we have a meaningful and concise as possible, but unique, name for each RefSeq.
    with open(args.output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(["accession", "tree_name"])
        for org, assemblies in organism_assemblies.items():
            org_san = sanitize_name(org)
            refseq_to_candidate = {}
            candidate_to_refseq = {}
            multiple_assemblies = (len(assemblies.keys()) > 1)
            dups = []
            for asm_acc, asm_info in assemblies.items():
                candidate_base = org_san
                if multiple_assemblies:
                    candidate_base += "_" + choose_assembly_suffix(asm_info)
                has_segments = len(asm_info.sequences) > 1
                for seq in asm_info.sequences:
                    candidate = candidate_base
                    if candidate_base.endswith("_"):
                        candidate += sanitize_name(seq.accession)
                        if has_segments and seq.segment:
                            candidate += "_" + sanitize_name(seq.segment)
                    elif has_segments:
                        suffix = sanitize_name(seq.segment if seq.segment else seq.accession)
                        candidate += "_" + suffix
                    # Check for duplicates; we'll fix those up later
                    if candidate in candidate_to_refseq:
                        if candidate_to_refseq[candidate] not in dups:
                            dups.append(candidate_to_refseq[candidate])
                        dups.append(seq.accession)
                    candidate_to_refseq[candidate] = seq.accession
                    refseq_to_candidate[seq.accession] = candidate
            for acc in dups:
                refseq_to_candidate[acc] += "_" + sanitize_name(acc)
            # By this point all should be resolved uniquely, add to refseq_to_name
            for acc, candidate in refseq_to_candidate.items():
                writer.writerow([acc, candidate])


if __name__ == "__main__":
    main()
