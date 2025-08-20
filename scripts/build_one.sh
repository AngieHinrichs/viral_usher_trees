#!/bin/bash

set -beEu -o pipefail

usage() {
    echo ""
    echo "usage: $0 trees_subdir"
    echo "trees_subdir is a subdirectory of trees/, e.g. measles.  Run this script in the repo top level."
    echo ""
}

if (( $# != 1 )); then
    usage
    exit 1
fi

subdir=$1

if [[ ! -d trees ]]; then
    echo ""
    echo "This script needs to be run in the top level of the repo, i.e. trees needs to be a subdirectory."
    echo ""
    exit 1
fi

if [[ ! -d trees/$subdir ]]; then
    echo ""
    echo "Can't find directory trees/$subdir"
    echo ""
    exit 1
fi

if [[ ! -s trees/$subdir/config.toml ]]; then
    echo ""
    echo "Can't find required viral_usher config file trees/$subdir/config.toml"
    echo ""
    exit 1
fi

# Config needs an absolute path; swap in 
sed -e "s@^workdir = '\.@workdir = '"$(pwd)"@" trees/$subdir/config.toml > trees/$subdir/local.toml
viral_usher build --config trees/$subdir/local.toml

rm -f trees/$subdir/{*.gzintermediate*,*.zip,*.fasta*,*.gbff,*.nh,empty*,*.log*,*.vcf.gz,data_report*,local.toml,changed_nodes,rename.tsv}
rm -f trees/$subdir/{nextclade.clade.tsv,mutation-paths.txt,placement_stats.tsv,usher_sampled.pb.gz,optimized.unfiltered.pb.gz}
