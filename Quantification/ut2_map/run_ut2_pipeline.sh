#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") <te1.nii|nii.gz> <te2.nii|nii.gz> <registration_ref|DICOM_dir> <segmentation_ref|DICOM_dir> <output_dir>

Example:
  $(basename "$0") te1.nii te2.nii /path/to/reg_ref.nii.gz /path/to/seg_ref.nii.gz out
USAGE
}

if [[ $# -ne 5 ]]; then
  usage >&2
  exit 1
fi

# Inputs passed to the pipeline.
TE1_INPUT=$1
TE2_INPUT=$2
REG_REF_INPUT=$3
SEG_REF_INPUT=$4
OUTDIR=$5

# Fixed rescaling settings.
RESCALING_LABELS="9"
RESCALING_THRESHOLD="170"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"

# Send all console output to both the terminal and a log file.
LOGFILE="$OUTDIR/pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

log "Starting UT2 pipeline"
log "Output directory: $OUTDIR"
log "Log file: $LOGFILE"
log "Step 1/3: registration"
"$SCRIPT_DIR/registration.sh" "$TE1_INPUT" "$TE2_INPUT" "$REG_REF_INPUT" "$OUTDIR"

log "Step 2/3: segmentation"
"$SCRIPT_DIR/segmentation.sh" "$SEG_REF_INPUT" "$RESCALING_LABELS" "$RESCALING_THRESHOLD" "$OUTDIR"

log "Step 3/3: UT2 subtraction"
"$SCRIPT_DIR/ut2_subtraction.sh" "$OUTDIR/te1_to_registration_ref.nii.gz" "$OUTDIR/te2_to_registration_ref_to_te1.nii.gz" "$OUTDIR/rescaling_mask_segmentation_ref.nii.gz" "$OUTDIR"

log "Pipeline completed successfully"
log "Final outputs are in: $OUTDIR"