#!/usr/bin/env bash
set -euo pipefail

# Build the rescaling mask used for UT2 subtraction by running MindGlide,
# selecting the requested labels, and applying the intensity threshold.
#
# Usage:
#   ./segmentation.sh segmentation_ref.nii.gz "1,2,3" 1200 path/to/outdir
#   ./segmentation.sh /path/to/dicom_dir "1 2 3" 1200 path/to/outdir

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.sh"
module load "$FSL_DIR"

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") <segmentation_ref|DICOM_dir> <rescaling_labels_csv> <rescaling_threshold> <output_dir>

Examples:
  $(basename "$0") segmentation_ref.nii.gz "1,2,3" 1200 out
  $(basename "$0") /path/to/dicom_dir "1 2 3" 1200 out
USAGE
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Error: required command not found: $1" >&2; exit 1; }
}

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

cleanup() {
  if [[ -n "${WORKDIR:-}" && -d "${WORKDIR:-}" ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

# Validate the required arguments and inputs.
if [[ $# -ne 4 ]]; then
  usage >&2
  exit 1
fi

SEG_REF_INPUT=$1
RESCALING_LABELS_RAW=$2
RESCALING_THRESHOLD=$3
OUTDIR=$4

[[ -e "$SEG_REF_INPUT" ]] || { echo "Error: segmentation reference not found: $SEG_REF_INPUT" >&2; exit 1; }

require_cmd fslreorient2std
require_cmd fslmaths
require_cmd dcm2niix
require_cmd bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
MINDGLIDE_SCRIPT="$SCRIPT_DIR/run_mindglide.sh"
[[ -x "$MINDGLIDE_SCRIPT" ]] || { echo "Error: expected executable mindglide script at $MINDGLIDE_SCRIPT" >&2; exit 1; }

# Create a temporary workspace for intermediate files.
mkdir -p "$OUTDIR"
WORKDIR=$(mktemp -d "$OUTDIR/.segmentation_tmp.XXXXXX")

log "Starting segmentation"
log "Segmentation reference: $SEG_REF_INPUT"
log "Rescaling labels: $RESCALING_LABELS_RAW"
log "Rescaling threshold: $RESCALING_THRESHOLD"

# Convert the segmentation reference to NIfTI if it is a DICOM directory.
if [[ -d "$SEG_REF_INPUT" ]]; then
  log "Converting DICOM segmentation reference with dcm2niix"
  dcm2niix -z y -o "$WORKDIR" -f segmentation_ref "$SEG_REF_INPUT" >/dev/null
  SEG_REF_NIFTI="$WORKDIR/segmentation_ref.nii.gz"
  [[ -f "$SEG_REF_NIFTI" ]] || { echo "Error: dcm2niix did not create $SEG_REF_NIFTI" >&2; exit 1; }
  cp -f "$SEG_REF_NIFTI" "$OUTDIR/segmentation_ref.nii.gz"
else
  SEG_REF_NIFTI="$SEG_REF_INPUT"
  cp -f "$SEG_REF_NIFTI" "$OUTDIR/segmentation_ref.nii.gz"
fi

# Reorient once before running mindglide and thresholding.
SEG_REF_REORIENT="$WORKDIR/segmentation_ref_reoriented.nii.gz"
log "Reorienting segmentation reference"
fslreorient2std "$SEG_REF_NIFTI" "$SEG_REF_REORIENT"

# Run your external mindglide wrapper script.
MINDGLIDE_OUT="$OUTDIR/mindglide_mask_segmentation_ref.nii.gz"
log "Running mindglide"
"$MINDGLIDE_SCRIPT" "$SEG_REF_REORIENT" "$MINDGLIDE_OUT"
[[ -f "$MINDGLIDE_OUT" ]] || { echo "Error: mindglide output not found: $MINDGLIDE_OUT" >&2; exit 1; }

# Parse labels from a comma- or whitespace-separated string.
LABELS_CSV=${RESCALING_LABELS_RAW// /,}
IFS=',' read -r -a LABELS <<< "$LABELS_CSV"
FILTERED_LABELS=()
for label in "${LABELS[@]}"; do
  [[ -n "$label" ]] && FILTERED_LABELS+=("$label")
done
[[ ${#FILTERED_LABELS[@]} -gt 0 ]] || { echo "Error: no rescaling labels provided" >&2; exit 1; }

# Build a binary mask that includes only the requested labels.
LABEL_MASK="$WORKDIR/rescaling_mask_labels.nii.gz"
FIRST_LABEL=1
for label in "${FILTERED_LABELS[@]}"; do
  log "Extracting label $label"
  LABEL_BIN="$WORKDIR/label_${label}.nii.gz"
  fslmaths "$MINDGLIDE_OUT" -thr "$label" -uthr "$label" -bin "$LABEL_BIN"
  if [[ $FIRST_LABEL -eq 1 ]]; then
    cp -f "$LABEL_BIN" "$LABEL_MASK"
    FIRST_LABEL=0
  else
    fslmaths "$LABEL_MASK" -add "$LABEL_BIN" -bin "$LABEL_MASK.tmp.nii.gz"
    mv -f "$LABEL_MASK.tmp.nii.gz" "$LABEL_MASK"
  fi
done

# Remove voxels in the segmentation reference above the threshold.
THRESH_MASK="$WORKDIR/below_threshold_mask.nii.gz"
log "Applying threshold filter"
fslmaths "$SEG_REF_REORIENT" -uthr "$RESCALING_THRESHOLD" -bin "$THRESH_MASK"

# Final rescaling mask is the label mask restricted by the threshold mask.
RESCALING_MASK="$OUTDIR/rescaling_mask_segmentation_ref.nii.gz"
log "Creating final rescaling mask"
fslmaths "$LABEL_MASK" -mas "$THRESH_MASK" "$RESCALING_MASK"

log "Segmentation step completed"
log "Outputs written to: $OUTDIR"
