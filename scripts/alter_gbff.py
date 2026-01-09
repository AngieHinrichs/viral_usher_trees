#!/usr/bin/env python3
"""
Given a GenBank flat file (GBFF) with viral genome sequences, the accession of one sequence
in the GBFF file, and a FASTA file with one sequence that has the same coordinates as the GBFF
sequence but some substituted bases, create a new GBFF file with the altered sequence and
annotations.
"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import sys
from typing import List


def read_fasta_one_sequence(fasta_path: str) -> SeqRecord:
    """Read a FASTA file and return the single SeqRecord it contains."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"FASTA file {fasta_path} must contain exactly one sequence.")
    return records[0]


def read_gbff_accession(gbff_path: str, accession: str) -> SeqRecord:
    """Read a GBFF file and return the SeqRecord with the specified accession."""
    for record in SeqIO.parse(gbff_path, "genbank"):
        if record.id == accession:
            return record
    raise ValueError(f"Accession {accession} not found in GBFF file {gbff_path}.")


def alter_gbff(gbff_record: SeqRecord, fasta_record: SeqRecord) -> SeqRecord:
    """Alter the sequence of the specified accession in the GBFF file."""
    # Create a new SeqRecord with the altered sequence
    altered_record = SeqRecord(
        seq=fasta_record.seq,
        id=fasta_record.name,
        name=fasta_record.name,
        description=fasta_record.description,
        annotations=gbff_record.annotations,
        features=gbff_record.features
    )
    for feature in altered_record.features:
        # If feature includes sequence, replace the sequence using the altered sequence
        if feature.type in {"CDS", "gene", "mRNA"}:
            if isinstance(feature.location, FeatureLocation):
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                feature_seq = altered_record.seq[start:end]
                if strand == -1:
                    feature_seq = feature_seq.reverse_complement()
                if feature.type == "CDS":
                    # Update the translation in the feature's qualifiers
                    translation = str(feature_seq.translate(to_stop=False))
                    feature.qualifiers["translation"] = [translation]
    return altered_record


def write_gbff(records: List[SeqRecord], output_path: str):
    """Write the list of SeqRecords to a GBFF file."""
    with open(output_path, "w") as output_handle:
        SeqIO.write(records, output_handle, "genbank")


def alter_gbff_file(gbff_file: str, accession: str, fasta_file: str, output_file: str):
    fasta_record = read_fasta_one_sequence(fasta_file)
    gbff_record = read_gbff_accession(gbff_file, accession)
    if len(gbff_record.seq) != len(fasta_record.seq):
        print("Error: The length of the FASTA sequence must match the length of the GBFF sequence.", file=sys.stderr)
        sys.exit(1)
    altered_record = alter_gbff(gbff_record, fasta_record)
    write_gbff([altered_record], output_file)


def main():
    parser = argparse.ArgumentParser(
        description="Given a GenBank flat file (GBFF) with viral genome sequences, "
                    "the accession of one sequence in the GBFF file, and a FASTA file "
                    "with one sequence that has the same coordinates as the GBFF "
                    "sequence but some substituted bases, create a new GBFF file with "
                    "the altered sequence and annotations."
    )
    parser.add_argument("--gbff", required=True, help="Input GenBank flat file")
    parser.add_argument("--accession", required=True, help="Accession of sequence in GBFF to be altered")
    parser.add_argument("--fasta", required=True, help="FASTA file with one sequence to use for alteration")
    parser.add_argument("--output", required=True, help="Output GenBank flat file")
    args = parser.parse_args()
    alter_gbff_file(args.gbff, args.accession, args.fasta, args.output)


if __name__ == "__main__":
    main()
