#!/bin/bash

set -beEu -o pipefail

usage() {
    echo ""
    echo "usage: $0 job_index"
    echo ""
    echo "The file job_list.<job_index>.txt must contain one tree subdir name per line."
    echo "This script must be run in the top-level of the repo."
    echo "It runs scripts/build_one.sh for each tree."
    echo ""
}

if (( $# != 1 )); then
    usage
    exit 1
fi

job_index=$1

if [[ ! -d trees ]]; then
    echo ""
    echo "This script needs to be run in the top level of the repo, i.e. trees needs to be a subdirectory."
    echo ""
    exit 1
fi

job_list=job_list.$job_index.txt

if [[ ! -e $job_list ]]; then
    echo ""
    echo "Error: Required file $job_list does not exist."
    echo ""
    exit 1
fi

fail_file=failed_trees.$job_index.txt
pass_file=built_trees.$job_index.txt
cp /dev/null $fail_file
cp /dev/null $pass_file

mkdir -p logs

while read tree_name; do
    echo ""
    if [[ -d trees/$tree_name && -s trees/$tree_name/config.toml ]]; then
        echo "Building $tree_name at $(date)"
        if ./scripts/build_one.sh $tree_name >& logs/$tree_name.log.txt; then
            echo $tree_name >> $pass_file
            echo "$tree_name built successfully."
        else
            echo $tree_name >> $fail_file
            echo "$tree_name failed, added to $fail_file"
        fi
    else
        echo "$tree_name doesn't have config_file, it should not be in $job_list"
    fi
done < $job_list
echo ""
echo "Completed all jobs successfully at $(date)."
wc -l $pass_file $fail_file

