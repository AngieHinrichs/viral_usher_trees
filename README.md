# viral_usher_trees

This repository stores [UShER](https://github.com/yatisht/usher) phylogenetic trees
for viral RefSeqs that have at least 1000 associated INSDC/GenBank sequences,
limited to trees with at least 100 near-full-length alignable sequences.
The trees are uilt by [viral_usher](https://github.com/AngieHinrichs/viral_usher) and
automatically updated on a monthly basis.

All tree files are in the `trees` directory, with a subdirectory per pathogen.
Each subdirectory contains the following files:
- **tree.jsonl.gz**: Taxonium tree file, viewable using https://taxonium.org/
- **viz.nwk.gz**: tree in Newick format; tip names include isolate, accession, collection date
- **viz.pb.gz**: tree in UShER protobuf format; tip names include isolate, accession, collection date
- **optimized.pb.gz**: tree in UShER protobuf format; tip names are INSDC/GenBank accessions
- **output_stats.tsv**: summary of sequence counts (available, filtered, aligned, in final tree)
- **metadata.tsv.gz**: metadata file describing input sequences in TSV format

Each pathogen subdirectory of `trees` also has a file `config.toml` that contains configuration parameters for the viral_usher build.

## Browse Trees Interactively
**[Browse/search trees](https://AngieHinrichs.github.io/viral_usher_trees/)**
