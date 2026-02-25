#!/bin/bash

set -beEu -o pipefail

while read tree method; do
    if [[ -d trees/$tree && -s trees/$tree/config.toml ]]; then
        refAcc=$(grep refseq_acc trees/$tree/config.toml  | cut -d\' -f 2)
        echo $tree: $method $refAcc
        if [[ $method == "treetime" ]]; then
            python scripts/config_reroot.py -t $tree -f ./trees/$tree/treetime_rerooted_$refAcc.fa -g ./trees/$tree/treetime_rerooted_$refAcc.gbff
            cp -p trees/$tree/timetree_rerooted.pb.gz trees/$tree/optimized.pb.gz
            git add trees/$tree/treetime_rerooted_$refAcc.fa trees/$tree/treetime_rerooted_$refAcc.gbff
        elif [[ $method == "midpoint" ]]; then
            if [[ -s trees/$tree/midpoint_rerooted.fasta ]]; then
                mv trees/$tree/midpoint_rerooted.fasta trees/$tree/midpoint_rerooted_$refAcc.fa
                mv trees/$tree/midpoint_rerooted.gbff trees/$tree/midpoint_rerooted_$refAcc.gbff
            fi
            python scripts/config_reroot.py -t $tree -f ./trees/$tree/midpoint_rerooted_$refAcc.fa -g ./trees/$tree/midpoint_rerooted_$refAcc.gbff
            cp -p trees/$tree/midpoint_rerooted.pb.gz trees/$tree/optimized.pb.gz
            git add trees/$tree/midpoint_rerooted_$refAcc.fa trees/$tree/midpoint_rerooted_$refAcc.gbff
        elif [[ $method == "outgroup" ]]; then
            python scripts/config_reroot.py -t $tree -f ./trees/$tree/outgroup_rerooted_$refAcc.fa -g ./trees/$tree/outgroup_rerooted_$refAcc.gbff
            cp -p trees/$tree/rerooted_no_outgroup_optimized.pb.gz trees/$tree/optimized.pb.gz
            git add trees/$tree/outgroup_rerooted_$refAcc.fa trees/$tree/outgroup_rerooted_$refAcc.gbff
        else
            -e "\nError: unrecognized method '$method'"
            false
        fi
        bash scripts/build_one.sh $tree
    fi
done
