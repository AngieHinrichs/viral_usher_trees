#!/usr/bin/env python3
"""
Save the sequences from a GenBank flat file (GBFF) to a FASTA file.
"""

import argparse
from Bio import SeqIO
import sys


def gbff_to_fasta(gbff_path: str, fasta_path: str):
    """Convert a GBFF file to a FASTA file."""
    try:
        with open(gbff_path, "r") as gbff_handle, open(fasta_path, "w") as fasta_handle:
            records = SeqIO.parse(gbff_handle, "genbank")
            count = SeqIO.write(records, fasta_handle, "fasta")
            print(f"Converted {count} records from {gbff_path} to {fasta_path}.", file=sys.stderr)
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Convert GBFF to FASTA")
    parser.add_argument("gbff", help="Path to the input GBFF file")
    parser.add_argument("fasta", help="Path to the output FASTA file")
    args = parser.parse_args()

    gbff_to_fasta(args.gbff, args.fasta)

if __name__ == "__main__":
    main()
