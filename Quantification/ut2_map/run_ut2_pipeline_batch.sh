#!/usr/bin/env bash
set -euo pipefail

# This script runs run_ut2_pipeline.sh for multiple subjects
# using an input text file organized in blocks.
#
# INPUT FILE FORMAT:
# The input file should contain groups of 5 lines per subject:
#
#   1. Path to TE1 image (.nii or .nii.gz)
#   2. Path to TE2 image (.nii or .nii.gz)
#   3. Path to registration reference (NIfTI file or DICOM directory)
#   4. Path to segmentation reference (NIfTI file or DICOM directory)
#   5. Output directory
#
# Each subject block must contain exactly 5 non-empty lines.
#
# Blocks should be separated by a blank line.
#
# Lines starting with '#' are treated as comments and ignored.
#
# EXAMPLE INPUT FILE:
# # Subject 1
# /data/sub1/te1.nii
# /data/sub1/te2.nii
# /data/sub1/reg_ref
# /data/sub1/seg_ref
# /data/sub1/output
#
# # Subject 2
# /data/sub2/te1.nii
# /data/sub2/te2.nii
# /data/sub2/reg_ref
# /data/sub2/seg_ref
# /data/sub2/output
#
# USAGE:
#   ./run_ut2_pipeline_batch.sh input_blocks.txt

usage() {
  echo "Usage: $0 <input_blocks.txt>"
  exit 1
}

# Validate the single required argument.
[[ $# -ne 1 ]] && usage

INPUT_FILE=$1
[[ -f "$INPUT_FILE" ]] || { echo "Error: file not found: $INPUT_FILE"; exit 1; }

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIPELINE="$SCRIPT_DIR/run_ut2_pipeline.sh"

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

log "Starting batch pipeline using block format"

# Accumulate one subject at a time until a blank line ends the block.
block=()
line_num=0
subject_id=0
failures=0

while IFS= read -r line || [[ -n "$line" ]]; do
  line_num=$((line_num + 1))

  # Skip comment lines
  [[ "$line" =~ ^# ]] && continue

  # Blank line triggers execution of current block
  if [[ -z "$line" ]]; then
    if [[ ${#block[@]} -gt 0 ]]; then
      subject_id=$((subject_id + 1))

      if [[ ${#block[@]} -ne 5 ]]; then
        log "ERROR: subject $subject_id has ${#block[@]} lines (expected 5)"
        failures=$((failures + 1))
      else
        log "Running subject $subject_id"
        if ! "$PIPELINE" "${block[@]}"; then
          log "ERROR: subject $subject_id failed"
          failures=$((failures + 1))
        fi
      fi

      block=()
    fi
    continue
  fi

  block+=("$line")

done < "$INPUT_FILE"

# Handle final block if file does not end with blank line
if [[ ${#block[@]} -gt 0 ]]; then
  subject_id=$((subject_id + 1))

  if [[ ${#block[@]} -ne 5 ]]; then
    log "ERROR: subject $subject_id has ${#block[@]} lines (expected 5)"
    failures=$((failures + 1))
  else
    log "Running subject $subject_id"
    if ! "$PIPELINE" "${block[@]}"; then
      log "ERROR: subject $subject_id failed"
      failures=$((failures + 1))
    fi
  fi
fi

log "Batch complete. Failures: $failures"

[[ $failures -gt 0 ]] && exit 1