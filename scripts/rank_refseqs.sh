#!/bin/bash

set -beEu -o pipefail

scriptDir=$(dirname "${BASH_SOURCE[0]}")

email=$USER'@soe.ucsc.edu'

mkdir -p rank_refseqs

# Get metadata for all viral RefSeq genomes from NCBI Virus
ncbiVirusUrl='https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B%21tag%3DSeqType_s%7DSeqType_s%3A%28%22Nucleotide%22%29&fq=SourceDB_s%3A%28%22RefSeq%22%29&q=%2A%3A%2A&cmd=download&dlfmt=csv&fl=accession%3AAccVer_s%2Corganism%3ATaxName_s%2Cassembly%3ASetAcc_s%2Cisolate%3AIsolate_s%2Cstrain%3AStrain_s%2Cserotype%3ASerotype_s%2Csegment%3ASegment_s%2Clength%3ASLen_i%2Chost%3AHost_s&sort=id+asc&email='$email

curl -Ss "$ncbiVirusUrl" | csvtk csv2tab > rank_refseqs/refseq_metadata.tsv

# Query NCBI to find the species-level Taxonomy ID for each RefSeq (takes ~15 minutes)
$scriptDir/refseq_to_taxid.py -i rank_refseqs/refseq_metadata.tsv -e $email -o rank_refseqs/refseq_taxid.tsv

# Query NCBI to find out how many GenBank sequences there are for each Taxonomy ID (takes ~5 hours)
$scriptDir/refseq_taxid_to_genbank_count.py -i rank_refseqs/refseq_taxid.tsv -e $email -o rank_refseqs/refseq_taxid_genbank_count.tsv

# Join it all together
csvtk join -t rank_refseqs/refseq_metadata.tsv rank_refseqs/refseq_taxid_genbank_count.tsv > rank_refseqs/refseq_metadata_plus.tsv

# Find RefSeqs for species that have at least 1000 sequences in GenBank (including the RefSeqs themselves)
awk -v 'FS=\t' -v 'OFS=\t' '$2 != "" && $11 >= 1000' rank_refseqs/refseq_metadata_plus.tsv \
| csvtk sort -t -k11:nr -k2 -k3 \
    > rank_refseqs/refseq_metadata_ranked.tsv

# Determine a meaningful-as-possible unique name for each RefSeq that makes the cut
# (Next time around, pass in the existing name mapping as -i to preserve it and only
#  assign new names to newly included RefSeqs)
$scriptDir/refseq_metadata_to_tree_name.py \
    -m rank_refseqs/refseq_metadata_ranked.tsv \
    -o rank_refseqs/refseq_to_tree_name.tsv
