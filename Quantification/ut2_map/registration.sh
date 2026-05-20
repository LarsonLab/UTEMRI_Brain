#!/usr/bin/env bash
set -euo pipefail

# Register the TE1 and TE2 images to a shared reference, then align
# TE2 into TE1 space for downstream UT2 subtraction.
#
# Usage:
#   ./registration.sh te1.nii te2.nii /path/to/reg_ref.nii.gz path/to/outdir

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.sh"
module load "$FSL_DIR"
module load "$ANTS_DIR"

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") <te1.nii|nii.gz> <te2.nii|nii.gz> <registration_ref|DICOM_dir> <output_dir>
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

# Validate the required arguments and input files.
if [[ $# -ne 4 ]]; then
  usage >&2
  exit 1
fi

TE1_INPUT=$1
TE2_INPUT=$2
REG_REF_INPUT=$3
OUTDIR=$4

for f in "$TE1_INPUT" "$TE2_INPUT"; do
  [[ -e "$f" ]] || { echo "Error: input not found: $f" >&2; exit 1; }
done
[[ -e "$REG_REF_INPUT" ]] || { echo "Error: registration reference not found: $REG_REF_INPUT" >&2; exit 1; }

require_cmd fslreorient2std
require_cmd flirt
require_cmd dcm2niix
require_cmd N4BiasFieldCorrection

# Create a temporary workspace for intermediate files.
mkdir -p "$OUTDIR"
WORKDIR=$(mktemp -d "$OUTDIR/.registration_tmp.XXXXXX")

log "Starting registration"
log "TE1: $TE1_INPUT"
log "TE2: $TE2_INPUT"
log "Registration reference: $REG_REF_INPUT"

# Convert the registration reference to NIfTI if it is a DICOM directory.
if [[ -d "$REG_REF_INPUT" ]]; then
  log "Converting DICOM registration reference with dcm2niix"
  dcm2niix -z y -o "$WORKDIR" -f registration_ref "$REG_REF_INPUT" >/dev/null
  REG_REF_NIFTI="$WORKDIR/registration_ref.nii.gz"
  [[ -f "$REG_REF_NIFTI" ]] || { echo "Error: dcm2niix did not create $REG_REF_NIFTI" >&2; exit 1; }
  cp -f "$REG_REF_NIFTI" "$OUTDIR/registration_ref.nii.gz"
else
  REG_REF_NIFTI="$REG_REF_INPUT"
  cp -f "$REG_REF_NIFTI" "$OUTDIR/registration_ref.nii.gz"
fi

# Reorient all three images to standard orientation before bias correction and registration.
TE1_REORIENT="$WORKDIR/te1_reoriented.nii.gz"
TE2_REORIENT="$WORKDIR/te2_reoriented.nii.gz"
REG_REF_REORIENT="$WORKDIR/registration_ref_reoriented.nii.gz"

log "Reorienting TE1, TE2, and registration reference"
fslreorient2std "$TE1_INPUT" "$TE1_REORIENT"
fslreorient2std "$TE2_INPUT" "$TE2_REORIENT"
fslreorient2std "$REG_REF_NIFTI" "$REG_REF_REORIENT"

run_n4_bias_correction() {
  local in_img=$1
  local out_prefix=$2
  local corrected_img="$WORKDIR/${out_prefix}_corrected.nii.gz"
  local biasfield_img="$WORKDIR/${out_prefix}_biasfield.nii.gz"

  # N4 is ANTs' standard bias-field correction method. The command below
  # keeps the behavior simple and robust, while still writing both the
  # corrected image and the estimated bias field for downstream inspection.
  N4BiasFieldCorrection \
    -d 3 \
    -i "$in_img" \
    -o "[$corrected_img,$biasfield_img]" \
    -s 2 \
    -b [200] \
    -c [50x50x30x20,0.0001] \
    -v

  [[ -f "$corrected_img" ]] || { echo "Error: N4 did not create $corrected_img" >&2; exit 1; }
  [[ -f "$biasfield_img" ]] || { echo "Error: N4 did not create $biasfield_img" >&2; exit 1; }

  mv -f "$corrected_img" "$OUTDIR/${out_prefix}.nii.gz"
  mv -f "$biasfield_img" "$OUTDIR/${out_prefix}map.nii.gz"
}

log "Bias-field correcting TE1 with ANTs N4"
run_n4_bias_correction "$TE1_REORIENT" "te1_bfc"

log "Bias-field correcting TE2 with ANTs N4"
run_n4_bias_correction "$TE2_REORIENT" "te2_bfc"

# Register each TE image to the reference, then register TE2 into TE1 space.
log "Registering TE1 to registration reference"
flirt -in "$OUTDIR/te1_bfc.nii.gz" -ref "$REG_REF_REORIENT" -out "$OUTDIR/te1_to_registration_ref.nii.gz"

log "Registering TE2 to registration reference"
flirt -in "$OUTDIR/te2_bfc.nii.gz" -ref "$REG_REF_REORIENT" -out "$OUTDIR/te2_to_registration_ref.nii.gz"

log "Registering TE2-in-reference-space to TE1-in-reference-space"
flirt -in "$OUTDIR/te2_to_registration_ref.nii.gz" -ref "$OUTDIR/te1_to_registration_ref.nii.gz" -out "$OUTDIR/te2_to_registration_ref_to_te1.nii.gz"

log "Registration step completed"
log "Outputs written to: $OUTDIR"
