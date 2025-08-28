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
        if timeout 120m ./scripts/build_one.sh $tree_name >& logs/$tree_name.log.txt; then
            echo $tree_name >> $pass_file
            echo "$tree_name built successfully."
        else
            echo $tree_name >> $fail_file
            echo "$tree_name failed, added to $fail_file"
            # In the context of a long run on a GitHub runner with little disk space,
            # remove all of the files generated for failed runs too, because they won't
            # be available after this job exits anyway.
            rm -f trees/$tree_name/{*.gzintermediate*,*.zip,*.fasta*,*.gbff,*.nh,empty*,*.log*,*.vcf.gz,data_report*,local.toml,changed_nodes,rename.tsv}
            rm -f trees/$tree_name/{nextclade.clade.tsv,mutation-paths.txt,placement_stats.tsv,usher_sampled.pb.gz,optimized.unfiltered.pb.gz}
            rm -f trees/$tree_name/{ncbi_virus_metadata.csv,tree_samples.txt}
            # Discard any changes to files under revision control
            changed_files=$(git status -s trees/$tree_name/ | awk '$1 == "M" {print $2;}')
            if [[ "$changed_files" != "" ]]; then
                echo "Discarding changes to files:"
                echo "$changed_files"
                git checkout HEAD -- $changed_files
            fi
            # Remove any orphaned new result files to prevent them from being added & committed later
            stray_files=$(git status -s trees/$tree_name/ | awk '$1 == "??" {print $2;}')
            if [[ "$stray_files" != "" ]]; then
                echo "Removing new files:"
                echo "$stray_files"
                git restore --staged $stray_files
                rm -f $stray_files
            fi
            # Add some delay in case the failure was a rejected connection or network failure
            sleep 30
        fi
    else
        echo "$tree_name doesn't have config_file, it should not be in $job_list"
    fi
done < $job_list
echo ""
echo "Completed all jobs successfully at $(date)."
wc -l $pass_file $fail_file

