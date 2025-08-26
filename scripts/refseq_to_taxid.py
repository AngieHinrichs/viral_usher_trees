#!/usr/bin/env python3
"""
Read a TSV file with RefSeq accessions and output a TSV mapping RefSeq accessions to taxonomy IDs
Initial version written by Claude, with a few tweaks (exit on query error, Unix newlines in output,
write output every batch in case of failure, query Taxonomy too to get species-level ID).
"""

import csv
import sys
import time
from pathlib import Path
from Bio import Entrez
import argparse
from typing import Dict, List

# Number of accessions per request to NCBI
default_batch_size = 100
# Delay between requests
default_delay = 1.0


def setup_entrez(email: str) -> None:
    """Setup Entrez with user email (required by NCBI)"""
    Entrez.email = email


def get_records(db, id_list):
    """Query NCBI and return parsed records.  May throw exceptions."""
    handle = Entrez.efetch(db=db, id=id_list, rettype="xml")
    records = Entrez.read(handle)
    handle.close()
    return records


def get_taxonomy_id_batch(accessions: List[str]) -> Dict[str, str]:
    """
    Get species-level taxonomy IDs for a batch of RefSeq accessions
    Returns a dictionary mapping accession -> taxonomy_id
    """
    if not accessions:
        return {}

    refseq_to_taxid = {}
    taxid_to_species_taxid = {}

    # Join accessions for batch query
    id_list = ",".join(accessions)

    try:
        # First, get the nucleotide records for these accessions and extract taxid
        records = get_records("nucleotide", id_list)
        for record in records:
            # Extract accession and taxonomy ID from the record
            accession = record.get('GBSeq_accession-version', record.get('GBSeq_locus', ''))
            # Look for taxonomy ID in organism info
            taxid = None
            for feature in record.get('GBSeq_feature-table', []):
                if feature.get('GBFeature_key') == 'source':
                    for qualifier in feature.get('GBFeature_quals', []):
                        if qualifier.get('GBQualifier_name') == 'db_xref':
                            value = qualifier.get('GBQualifier_value', '')
                            if value.startswith('taxon:'):
                                taxid = value.split(':')[1]
                                break
                    if taxid:
                        break
            if accession and taxid:
                refseq_to_taxid[accession] = taxid
                taxid_to_species_taxid[taxid] = taxid
            else:
                refseq_to_taxid[accession] = "N/A"

        # Next, query Taxonomy to get the species-level taxonomy id for each RefSeq-associated id
        id_list = ",".join(taxid_to_species_taxid.keys())
        records = get_records("taxonomy", id_list)
        for record in records:
            taxid = record.get("TaxId")
            rank = record.get("Rank", "no rank")
            lineage = record.get("LineageEx", [])
            if taxid and rank and lineage:
                # If taxid's rank is not species then trace back through lineage to find species-level taxid
                if rank != "species":
                    for ancestor in reversed(lineage):
                        if ancestor["Rank"] == "species":
                            taxid_to_species_taxid[taxid] = ancestor["TaxId"]
                            break

        # Finally, substitute all of the RefSeq-associated taxids with their species-level taxids
        for refseq, taxid in refseq_to_taxid.items():
            species_level_taxid = taxid_to_species_taxid.get(taxid)
            if species_level_taxid:
                refseq_to_taxid[refseq] = species_level_taxid

    except Exception as e:
        print(f"Error querying batch: {e}")
        sys.exit(1)

    return refseq_to_taxid


def process_tsv_file(input_file: str, output_file: str, email: str,
                     accession_column: str = None, batch_size: int = default_batch_size,
                     delay: float = default_delay) -> None:
    """
    Process TSV file and map RefSeq accessions to taxonomy IDs
    Args:
        input_file: Path to input TSV file
        output_file: Path to output TSV file
        email: Email for NCBI Entrez (required)
        accession_column: Name of column containing accessions (auto-detect if None)
        batch_size: Number of accessions to query at once
        delay: Delay between batches in seconds
    """

    setup_entrez(email)

    # Read input file
    print(f"Reading input file: {input_file}")
    with open(input_file, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)

    print(f"Found {len(rows)} rows with columns: {fieldnames}")

    # Auto-detect accession column if not specified
    if accession_column is None:
        # Look for columns that might contain RefSeq accessions
        possible_columns = [col for col in fieldnames if any(keyword in col.lower()
                            for keyword in ['accession', 'refseq', 'nc_', 'id'])]

        if possible_columns:
            accession_column = possible_columns[0]
            print(f"Auto-detected accession column: '{accession_column}'")
        else:
            # Default to first column
            accession_column = fieldnames[0]
            print(f"Using first column as accession column: '{accession_column}'")

    # Extract accessions
    accessions = []
    for row in rows:
        acc = row.get(accession_column, '').strip()
        if acc and (acc.startswith('NC_') or acc.startswith('AC_')):
            accessions.append(acc)

    print(f"Found {len(accessions)} RefSeq accessions to process")

    if not accessions:
        print("No RefSeq accessions found! Check your input file and column specification.\n" +
              "Use -c/--column to specify accession column.")
        return

    # Process in batches, writing output file as we go so if batch fails we get partial results
    total_batches = (len(accessions) + batch_size - 1) // batch_size

    print(f"\nWriting results to: {output_file}")
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        # Write header
        writer.writerow(['RefSeq_Accession', 'Taxonomy_ID'])
        successful = 0
        for i in range(0, len(accessions), batch_size):
            batch_num = (i // batch_size) + 1
            batch = accessions[i:i+batch_size]

            print(f"\nProcessing batch {batch_num}/{total_batches} ({len(batch)} accessions)")
            batch_results = get_taxonomy_id_batch(batch)
            for accession in batch:
                taxid = batch_results.get(accession, 'NOT_FOUND')
                writer.writerow([accession, taxid])
                if taxid not in ['N/A', 'NOT_FOUND']:
                    successful += 1

            # Add delay between batches to be nice to NCBI servers
            if i + batch_size < len(accessions):
                print(f"Waiting {delay} seconds before next batch...")
                time.sleep(delay)

    print("\nProcessing complete!")
    print(f"Successfully mapped {successful}/{len(accessions)} accessions")
    print(f"Results written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Map RefSeq accessions to NCBI Taxonomy IDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python refseq_to_taxid.py -i input.tsv -o output.tsv -e your.email@example.com
  python refseq_to_taxid.py -i data.tsv -o results.tsv -e user@domain.com -c "Accession" -b 10 -d 2.0
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file containing RefSeq accessions')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file for accession -> taxonomy ID mapping')
    parser.add_argument('-e', '--email', required=True,
                        help='Your email address (required by NCBI)')
    parser.add_argument('-c', '--column',
                        help='Column name containing RefSeq accessions (auto-detect if not specified)')
    parser.add_argument('-b', '--batch-size', type=int, default=default_batch_size,
                        help=f"Number of accessions to query per batch (default: {default_batch_size})")
    parser.add_argument('-d', '--delay', type=float, default=default_delay,
                        help=f"Delay between batches in seconds (default: {default_delay})")

    args = parser.parse_args()

    # Validate input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found")
        sys.exit(1)

    # Process the file
    try:
        process_tsv_file(
            input_file=args.input,
            output_file=args.output,
            email=args.email,
            accession_column=args.column,
            batch_size=args.batch_size,
            delay=args.delay
        )
    except KeyboardInterrupt:
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
