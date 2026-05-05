#!/usr/bin/env bash
set -euo pipefail

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

if [[ $# -ne 4 ]]; then
  usage >&2
  exit 1
fi

# Load FSL before checking for FSL tools.
FSL_DIR="SCS/fsl/fsl_latest"
module load "$FSL_DIR"

TE1_INPUT=$1
TE2_INPUT=$2
REG_REF_INPUT=$3
OUTDIR=$4

for f in "$TE1_INPUT" "$TE2_INPUT"; do
  [[ -e "$f" ]] || { echo "Error: input not found: $f" >&2; exit 1; }
done
[[ -e "$REG_REF_INPUT" ]] || { echo "Error: registration reference not found: $REG_REF_INPUT" >&2; exit 1; }

require_cmd fslreorient2std
require_cmd fast
require_cmd flirt
require_cmd dcm2niix

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

# FAST performs bias-field correction and also writes bias-field outputs.
log "Bias-field correcting TE1"
fast -B -b -o "$WORKDIR/te1_bfc" "$TE1_REORIENT"
log "Bias-field correcting TE2"
fast -B -b -o "$WORKDIR/te2_bfc" "$TE2_REORIENT"

# FAST output names can vary by version, so normalize them to the requested names.
for base in te1_bfc te2_bfc; do
  if [[ -f "$WORKDIR/${base}_restore.nii.gz" ]]; then
    mv -f "$WORKDIR/${base}_restore.nii.gz" "$OUTDIR/${base}.nii.gz"
  elif [[ -f "$WORKDIR/${base}.nii.gz" ]]; then
    mv -f "$WORKDIR/${base}.nii.gz" "$OUTDIR/${base}.nii.gz"
  else
    echo "Error: could not find restored output for $base" >&2
    exit 1
  fi

  if [[ -f "$WORKDIR/${base}_bias.nii.gz" ]]; then
    mv -f "$WORKDIR/${base}_bias.nii.gz" "$OUTDIR/${base}map.nii.gz"
  elif [[ -f "$WORKDIR/${base}_biasfield.nii.gz" ]]; then
    mv -f "$WORKDIR/${base}_biasfield.nii.gz" "$OUTDIR/${base}map.nii.gz"
  else
    echo "Error: could not find bias field output for $base" >&2
    exit 1
  fi

done

# Register each TE image to the reference, then register TE2 into TE1 space.
log "Registering TE1 to registration reference"
flirt -in "$OUTDIR/te1_bfc.nii.gz" -ref "$REG_REF_REORIENT" -out "$OUTDIR/te1_to_registration_ref.nii.gz"

log "Registering TE2 to registration reference"
flirt -in "$OUTDIR/te2_bfc.nii.gz" -ref "$REG_REF_REORIENT" -out "$OUTDIR/te2_to_registration_ref.nii.gz"

log "Registering TE2-in-reference-space to TE1-in-reference-space"
flirt -in "$OUTDIR/te2_to_registration_ref.nii.gz" -ref "$OUTDIR/te1_to_registration_ref.nii.gz" -out "$OUTDIR/te2_to_registration_ref_to_te1.nii.gz"

log "Registration step completed"
log "Outputs written to: $OUTDIR"
