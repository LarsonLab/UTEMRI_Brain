#!/usr/bin/env bash
# 
# Batch wrapper for lesion_segmentation.sh. This script reads a text file containing 
# one or more subject blocks and runs lesion_segmentation.sh once for each subject.
#
# INPUT FILE FORMAT
# Each subject is defined by 4 consecutive non-comment lines:
#   1. MPRAGE input: NIfTI (.nii or .nii.gz), or DICOM series directory
#   2. FLAIR input: NIfTI (.nii or .nii.gz), or DICOM series directory
#   3. Registration reference image: NIfTI (.nii or .nii.gz)
#   4. Output directory
#
# Subject blocks are separated by one or more blank lines.
# Lines beginning with '#' are treated as comments and ignored.
#
# EXAMPLE INPUT FILE
#   # Subject 001
#   /data/sub001/mprage
#   /data/sub001/flair
#   /refs/registration_ref.nii.gz
#   /results/sub001
#
#   # Subject 002
#   /data/sub002/mprage.nii.gz
#   /data/sub002/flair.nii.gz
#   /refs/registration_ref.nii.gz
#   /results/sub002
#
# USAGE
#   lesion_segmentation_batch.sh <subject_list.txt> <lesion_segmentation.sh>

set -euo pipefail

usage() {
    cat <<EOF
Usage:
    $(basename "$0") <subject_list.txt> <lesion_segmentation.sh>

Arguments:
    subject_list.txt       Subject input list
    lesion_segmentation.sh Single-subject processing script

EOF
}

if [[ $# -ne 2 ]]; then
    usage
    exit 1
fi

LIST="$1"
PIPELINE="$2"

if [[ ! -f "$LIST" ]]; then
    echo "ERROR: Subject list not found: $LIST" >&2
    exit 1
fi

if [[ ! -f "$PIPELINE" ]]; then
    echo "ERROR: Pipeline script not found: $PIPELINE" >&2
    exit 1
fi

block=()
subject_count=0

run_block() {
    [[ ${#block[@]} -eq 0 ]] && return

    if [[ ${#block[@]} -ne 4 ]]; then
        echo "ERROR: Expected 4 lines per subject block, got ${#block[@]}" >&2
        exit 1
    fi

    local mprage="${block[0]}"
    local flair="${block[1]}"
    local ref="${block[2]}"
    local outdir="${block[3]}"

    subject_count=$((subject_count + 1))

    echo
    echo "============================================================"
    echo "Subject ${subject_count}"
    echo "  MPRAGE : $mprage"
    echo "  FLAIR  : $flair"
    echo "  REF    : $ref"
    echo "  OUTDIR : $outdir"
    echo "============================================================"

    bash "$PIPELINE" \
        "$mprage" \
        "$flair" \
        "$ref" \
        "$outdir"

    block=()
}

while IFS= read -r line || [[ -n "$line" ]]; do

    # Remove Windows carriage returns if present
    line="${line%$'\r'}"

    # Skip comments
    [[ "$line" =~ ^[[:space:]]*# ]] && continue

    # Blank line marks end of a subject block
    if [[ -z "$line" ]]; then
        run_block
        continue
    fi

    block+=("$line")

done < "$LIST"

# Process final block
run_block

echo
echo "Finished processing ${subject_count} subject(s)."