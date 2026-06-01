#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") <mprage_in> <flair_in> <registration_ref.nii.gz> <output_dir>

Inputs:
  mprage_in            NIfTI file (.nii or .nii.gz) or DICOM series directory
  flair_in             NIfTI file (.nii or .nii.gz) or DICOM series directory
  registration_ref     NIfTI file (.nii or .nii.gz)
  output_dir           Directory where outputs will be saved

Outputs:
  mprage_to_registration_ref.nii.gz
  flair_to_mprage_to_registration_ref.nii.gz
  seg.nii.gz
  seg_to_registration_ref.nii.gz
  lesion_mask_to_registration_ref.nii.gz
EOF
}

module load "SCS/fsl/fsl_latest"
module load "SCS/freesurfer/7.4.1"

NCORES=${SLURM_CPUS_PER_TASK:-$(nproc)}

if [[ $# -ne 4 ]]; then
  usage
  exit 1
fi

MPRAGE_IN="$1"
FLAIR_IN="$2"
REF_IN="$3"
OUTDIR="$4"

mkdir -p "$OUTDIR"
WORKDIR="$OUTDIR/work"
mkdir -p "$WORKDIR"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Required command not found: $1" >&2
    exit 1
  }
}

for cmd in dcm2niix fslreorient2std flirt mri_convert run_samseg fslmaths; do
  require_cmd "$cmd"
done

if [[ ! -f "$REF_IN" ]]; then
  echo "ERROR: registration_ref must be a NIfTI file." >&2
  exit 1
fi

convert_input_to_nifti() {
  local inpath="$1"
  local tag="$2"
  local outpath="$WORKDIR/${tag}.nii.gz"

  if [[ -d "$inpath" ]]; then
    local dcm_out="$WORKDIR/${tag}_dcm2niix"
    mkdir -p "$dcm_out"

    # Convert the DICOM series directory to NIfTI.
    dcm2niix -z y -f "$tag" -o "$dcm_out" "$inpath" >/dev/null

    local converted
    converted="$(find "$dcm_out" -maxdepth 1 -type f \( -name "*.nii" -o -name "*.nii.gz" \) | head -n 1 || true)"

    if [[ -z "$converted" ]]; then
      echo "ERROR: dcm2niix did not produce a NIfTI for $inpath" >&2
      exit 1
    fi

    cp "$converted" "$outpath"
  else
    if [[ ! -f "$inpath" ]]; then
      echo "ERROR: Input not found: $inpath" >&2
      exit 1
    fi
    cp "$inpath" "$outpath"
  fi

  echo "$outpath"
}

reorient_std() {
  local inpath="$1"
  local tag="$2"
  local outpath="$WORKDIR/${tag}_std.nii.gz"
  fslreorient2std "$inpath" "$outpath"
  echo "$outpath"
}

echo "Preparing inputs..."
MPRAGE_NII="$(convert_input_to_nifti "$MPRAGE_IN" "mprage")"
FLAIR_NII="$(convert_input_to_nifti "$FLAIR_IN" "flair")"
REF_NII="$(convert_input_to_nifti "$REF_IN" "registration_ref")"

echo "Reorienting to standard space..."
MPRAGE_STD="$(reorient_std "$MPRAGE_NII" "mprage")"
FLAIR_STD="$(reorient_std "$FLAIR_NII" "flair")"
REF_STD="$(reorient_std "$REF_NII" "registration_ref")"

echo "Registering FLAIR to MPRAGE..."
FLAIR_TO_MPRAGE="$OUTDIR/flair_to_mprage_to_registration_ref_intermediate.nii.gz"
FLAIR_TO_MPRAGE_MAT="$OUTDIR/flair_to_mprage.mat"
flirt \
  -in "$FLAIR_STD" \
  -ref "$MPRAGE_STD" \
  -cost mutualinfo \
  -omat "$FLAIR_TO_MPRAGE_MAT" \
  -out "$FLAIR_TO_MPRAGE"

echo "Running SAMSEG..."
run_samseg \
  --input "$MPRAGE_STD" "$FLAIR_TO_MPRAGE" \
  --lesion \
  --lesion-mask-pattern 0 1 \
  --pallidum-separate \
  --threshold 0.3 \
  --output "$OUTDIR" \
  --threads "$NCORES"

SEG_MGZ="$OUTDIR/seg.mgz"
SEG_NII="$OUTDIR/seg.nii.gz"

if [[ ! -f "$SEG_MGZ" ]]; then
  echo "ERROR: SAMSEG output not found: $SEG_MGZ" >&2
  exit 1
fi

echo "Converting seg.mgz to NIfTI..."
mri_convert "$SEG_MGZ" "$SEG_NII"

echo "Registering MPRAGE to registration_ref..."
MPRAGE_TO_REF_MAT="$OUTDIR/mprage_to_registration_ref.mat"
MPRAGE_TO_REF="$OUTDIR/mprage_to_registration_ref.nii.gz"
flirt \
  -in "$MPRAGE_STD" \
  -ref "$REF_STD" \
  -cost mutualinfo \
  -omat "$MPRAGE_TO_REF_MAT" \
  -out "$MPRAGE_TO_REF"

echo "Applying MPRAGE->REF transform to flair and seg..."
FLAIR_TO_MPRAGE_TO_REF="$OUTDIR/flair_to_mprage_to_registration_ref.nii.gz"
SEG_TO_REF="$OUTDIR/seg_to_registration_ref.nii.gz"

flirt \
  -in "$FLAIR_TO_MPRAGE" \
  -ref "$REF_STD" \
  -applyxfm \
  -init "$MPRAGE_TO_REF_MAT" \
  -interp trilinear \
  -out "$FLAIR_TO_MPRAGE_TO_REF"

flirt \
  -in "$SEG_NII" \
  -ref "$REF_STD" \
  -applyxfm \
  -init "$MPRAGE_TO_REF_MAT" \
  -interp nearestneighbour \
  -out "$SEG_TO_REF"

echo "Extracting lesion label 99..."
LESION_MASK="$OUTDIR/lesion_mask_to_registration_ref.nii.gz"
fslmaths "$SEG_TO_REF" -thr 99 -uthr 99 -bin "$LESION_MASK"

echo "Done."
echo "Outputs:"
echo "  $MPRAGE_TO_REF"
echo "  $FLAIR_TO_MPRAGE_TO_REF"
echo "  $SEG_NII"
echo "  $SEG_TO_REF"
echo "  $LESION_MASK"