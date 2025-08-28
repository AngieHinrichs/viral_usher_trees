#!/usr/bin/env python3
"""
Read a TSV file with RefSeq accessions and Taxonomy IDs,
query NCBI for GenBank sequence counts, and output results to TSV.
Initial version written by Claude, with a few tweaks (arg formats, Unix newlines in output,
quintuple delays between retries, write output every batch in case of failure)
"""

import csv
import time
import sys
from typing import Tuple, List
from Bio import Entrez
import argparse
from pathlib import Path

# Delay between requests
default_delay = 1.0
# Default max number of tries if a query fails
default_max_tries = 5


def read_input_tsv(filepath: str) -> List[Tuple[str, str]]:
    """
    Read the input TSV file and return list of (RefSeq_accession, Taxonomy_ID) tuples.
    Args:
        filepath: Path to the input TSV file
    Returns:
        List of tuples containing (RefSeq accession, Taxonomy ID)
    """
    data = []
    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as file:
            # Try to detect if there's a header
            sample = file.read(1024)
            file.seek(0)
            sniffer = csv.Sniffer()
            has_header = sniffer.has_header(sample)

            reader = csv.reader(file, delimiter='\t')

            if has_header:
                next(reader)  # Skip header row

            for row in reader:
                if len(row) >= 2:
                    refseq_acc = row[0].strip()
                    taxonomy_id = row[1].strip()
                    data.append((refseq_acc, taxonomy_id))
    except FileNotFoundError:
        print(f"Error: Input file '{filepath}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    return data


def get_genbank_count(taxonomy_id: str, max_retries: int = default_max_tries) -> int:
    """
    Query NCBI to get the number of GenBank sequences for a given Taxonomy ID.
    Args:
        taxonomy_id: NCBI Taxonomy ID
        max_retries: Maximum number of retry attempts
    Returns:
        Number of GenBank sequences for the taxonomy ID
    """
    delay = default_delay
    for attempt in range(max_retries):
        try:
            # Search GenBank database for sequences with the given taxonomy ID
            search_term = f"txid{taxonomy_id}[Organism]"
            handle = Entrez.esearch(
                db="nucleotide",  # GenBank nucleotide database
                term=search_term,
                retmax=0  # We only want the count, not the actual records
            )
            search_results = Entrez.read(handle)
            handle.close()
            count = int(search_results["Count"])
            return count
        except Exception as e:
            print(f"Attempt {attempt + 1} failed for taxonomy ID {taxonomy_id}: {e}")
            if attempt < max_retries - 1:
                # Wait before retrying, quintuple delay between retries
                time.sleep(delay)
                delay *= 5
            else:
                print(f"Failed to get count for taxonomy ID {taxonomy_id} after {max_retries} attempts")
                return -1  # Return -1 to indicate failure
    return -1


def main():
    parser = argparse.ArgumentParser(description="Query NCBI for GenBank sequence counts by Taxonomy ID")
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file with RefSeq accessions and Taxonomy IDs')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file with RefSeq accessions, Taxonomy IDs and Genbank counts')
    parser.add_argument('-e', '--email', required=True,
                        help='Your email address (required by NCBI)')
    parser.add_argument('-d', '--delay', type=float, default=default_delay,
                        help=f"Delay between queries in seconds (default: {default_delay})")
    args = parser.parse_args()

    # Set the email for NCBI Entrez
    Entrez.email = args.email

    # Validate input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found")
        sys.exit(1)

    print(f"Reading input file: {args.input}")
    input_data = read_input_tsv(args.input)

    if not input_data:
        print("No data found in input file.")
        sys.exit(1)

    print(f"Found {len(input_data)} records to process")
    print(f"Writing output to {args.output}")
    with open(args.output, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n')
        writer.writerow(['RefSeq_Accession', 'Taxonomy_ID', 'GenBank_Count'])
        # Process each taxonomy ID
        unique_taxonomy_ids = {}  # Cache to avoid duplicate queries

        for i, (refseq_acc, taxonomy_id) in enumerate(input_data, 1):
            print(f"Processing {i}/{len(input_data)}: {refseq_acc} (Taxonomy ID: {taxonomy_id})")

            # Check if we've already queried this taxonomy ID
            if taxonomy_id in unique_taxonomy_ids:
                genbank_count = unique_taxonomy_ids[taxonomy_id]
            else:
                genbank_count = get_genbank_count(taxonomy_id)
                unique_taxonomy_ids[taxonomy_id] = genbank_count
                if genbank_count >= 0:
                    print(f"  Found {genbank_count} GenBank sequences")
                else:
                    print(f"  Failed to retrieve count for {taxonomy_id}")
                    sys.exit(1)
                # Add delay to be respectful to NCBI servers
                time.sleep(args.delay)
            writer.writerow([refseq_acc, taxonomy_id, genbank_count])

    # Print summary
    unique_taxa = len(unique_taxonomy_ids)

    print("\nSummary:")
    print(f"  Total records processed: {len(input_data)}")
    print(f"  Unique taxonomy IDs: {unique_taxa}")


if __name__ == "__main__":
    main()
