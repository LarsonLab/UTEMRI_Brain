#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") <te1.nii|nii.gz> <te2.nii|nii.gz> <registration_ref|DICOM_dir> <segmentation_ref|DICOM_dir> <output_dir>

Example:
  $(basename "$0") te1.nii te2.nii /path/to/reg_ref.nii.gz /path/to/seg_ref.nii.gz out "1,2,3" 1200
USAGE
}

if [[ $# -ne 7 ]]; then
  usage >&2
  exit 1
fi

TE1_INPUT=$1
TE2_INPUT=$2
REG_REF_INPUT=$3
SEG_REF_INPUT=$4
OUTDIR=$5

RESCALING_LABELS="9"
RESCALING_THRESHOLD="170"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "$OUTDIR"

"$SCRIPT_DIR/registration.sh" "$TE1_INPUT" "$TE2_INPUT" "$REG_REF_INPUT" "$OUTDIR"
"$SCRIPT_DIR/segmentation.sh" "$SEG_REF_INPUT" "$RESCALING_LABELS" "$RESCALING_THRESHOLD" "$OUTDIR"
"$SCRIPT_DIR/ut2_subtraction.sh" "$OUTDIR/te1_to_registration_ref.nii.gz" "$OUTDIR/te2_to_registration_ref_to_te1.nii.gz" "$OUTDIR/rescaling_mask_segmentation_ref.nii.gz" "$OUTDIR"

echo "Pipeline completed. Final outputs are in: $OUTDIR"
