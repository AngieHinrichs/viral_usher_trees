# viral_usher_trees

This repository stores [UShER](https://github.com/yatisht/usher) phylogenetic trees for common viral pathogens,
built by [viral_usher](https://github.com/AngieHinrichs/viral_usher) and
automatically updated on a monthly basis.

All tree files are in the `trees` directory, with a subdirectory per pathogen.
Each subdirectory contains the following files:
- optimized.pb.gz: tree in UShER protobuf format
- metadata.tsv.gz: metadata file describing input sequences in TSV format
- tree.jsonl.gz: Taxonium tree file, viewable using https://taxonium.org/
