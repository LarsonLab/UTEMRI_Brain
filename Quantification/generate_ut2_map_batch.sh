#!/bin/bash
set -euo pipefail

# Script to batch run generate_ut2_map.sh
# Usage: ./batch_run.sh job_list.txt path/to/main_script.sh
# job_list.txt should have one job per line, formatted as:
# <mprage_folder> <te1_file> <te2_file> <output_dir>

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 job_list.txt path/to/main_script.sh"
    exit 1
fi

JOB_LIST="$1"
MAIN_SCRIPT="$2"

if [[ ! -f "$JOB_LIST" ]]; then
    echo "Error: job list $JOB_LIST not found."
    exit 1
fi

if [[ ! -x "$MAIN_SCRIPT" ]]; then
    echo "Error: main script $MAIN_SCRIPT is not executable."
    echo "Try: chmod +x $MAIN_SCRIPT"
    exit 1
fi

# Loop through each line in job list
while read -r MPRAGE TE1 TE2 OUTDIR; do
    # Skip empty lines or comments
    [[ -z "$MPRAGE" || "$MPRAGE" =~ ^# ]] && continue

    echo "---------------------------------------"
    echo "Running job with:"
    echo "MPRAGE: $MPRAGE"
    echo "TE1:    $TE1"
    echo "TE2:    $TE2"
    echo "OUTDIR: $OUTDIR"
    echo "---------------------------------------"

    "$MAIN_SCRIPT" -m "$MPRAGE" -te1 "$TE1" -te2 "$TE2" -o "$OUTDIR"
done < "$JOB_LIST"

echo "All jobs finished."
