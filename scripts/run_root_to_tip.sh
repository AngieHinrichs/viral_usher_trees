#!/bin/bash
set -eu

for d in trees/*; do
    if [[ -d $d ]]; then
        v=$(basename $d)
        pushd $d
        if [[ -e timetree_rerooted.pb.gz ]]; then
            maybe_b="-b timetree_rerooted.pb.gz"
        else
            maybe_b=""
        fi
        if [[ -e rerooted_no_outgroup_optimized.pb.gz ]]; then
            maybe_c="-c rerooted_no_outgroup_optimized.pb.gz"
        else
            maybe_c=""
        fi
        python ../../scripts/root_to_tip.py \
               -m metadata.tsv.gz \
               -a midpoint_rerooted.pb.gz \
               -A accession \
               --title_a "midpoint" \
               $maybe_b \
               -B accession \
               --title_b "treetime" \
               $maybe_c \
               -C accession \
               --title_c "outgroup" \
               -o all.rtt \
               -t $v
        popd
    fi
done
