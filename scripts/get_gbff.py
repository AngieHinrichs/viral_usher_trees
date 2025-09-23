#!/usr/bin/env python3
"""
Given an accession, fetch the corresponding GenBank flat file (GBFF) from NCBI.
"""
import argparse
from Bio import Entrez, SeqIO
import sys


def fetch_gbff(accession: str, email: str) -> SeqIO.SeqRecord:
    """Fetch a GenBank flat file (GBFF) from NCBI given an accession."""
    Entrez.email = email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")


def main():
    parser = argparse.ArgumentParser(description="Fetch a GenBank flat file (GBFF) from NCBI given an accession.")
    parser.add_argument("accession", help="The accession of the sequence to fetch.")
    parser.add_argument("email", help="Email address to use with NCBI Entrez.")
    args = parser.parse_args()

    try:
        gbff_record = fetch_gbff(args.accession, args.email)
        with open(args.accession + ".gbff", "w") as out_handle:
            SeqIO.write(gbff_record, out_handle, "genbank")
    except Exception as e:
        print(f"Error fetching accession {args.accession}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
