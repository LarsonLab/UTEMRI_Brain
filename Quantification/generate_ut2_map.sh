#!/bin/bash
set -euo pipefail

# Constants
TE2_NORMALIZATION_FACTOR=1.38  # scaling factor to normalize TE2 image

# Load FSL module
FSL_DIR="SCS/fsl/fsl_latest"  # set to your FSL module path
module load "$FSL_DIR"

# Initialize variables
MPRAGE_DICOM_FOLDER=""
TE1_IN=""
TE2_IN=""
OUTDIR=""

# Parse command line flags
while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--mprage) MPRAGE_DICOM_FOLDER="$2"; shift 2 ;;
    -te1|--te1) TE1_IN="$2"; shift 2 ;;
    -te2|--te2) TE2_IN="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -m <MPRAGE_dicom_folder> -te1 <te1.nii.gz> -te2 <te2.nii.gz> -o <output_directory>"
      exit 0 ;;
    *) echo "Unknown option: $1"; echo "Use -h for help."; exit 1 ;;
  esac
done

# Validate inputs
if [[ -z "$MPRAGE_DICOM_FOLDER" || -z "$TE1_IN" || -z "$TE2_IN" || -z "$OUTDIR" ]]; then
  echo "ERROR: Missing required arguments."
  echo "Usage: $0 -m <MPRAGE_dicom_folder> -te1 <te1.nii.gz> -te2 <te2.nii.gz> -o <output_directory>"
  exit 1
fi

# Setup directories
mkdir -p "${OUTDIR}"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Working directory: ${WORKDIR}"
echo "Output directory: ${OUTDIR}"

# File naming conventions
MPRAGE_BASE="${OUTDIR}/MPRAGE"
TE1_BFC="${WORKDIR}/te1_bfc.nii.gz"
TE2_BFC="${WORKDIR}/te2_bfc.nii.gz"

# Convert MPRAGE DICOMs to gzipped NIfTI
echo "Converting MPRAGE DICOMs to NIfTI (.nii.gz)..."
dcm2niix -z y -f MPRAGE -o "${WORKDIR}" "${MPRAGE_DICOM_FOLDER}"

MPRAGE_NII="$(ls -1 "${WORKDIR}"/*.nii.gz 2>/dev/null | head -n1 || true)"
if [ -z "${MPRAGE_NII}" ]; then
  echo "ERROR: dcm2niix did not produce a .nii.gz file in ${WORKDIR}."
  exit 1
fi
cp "${MPRAGE_NII}" "${MPRAGE_BASE}.nii.gz"
MPRAGE_NII="${MPRAGE_BASE}.nii.gz"

# Reorient all inputs to standard orientation
echo "Reorienting all input images to standard orientation..."
REORIENTED_MPRAGE="${WORKDIR}/MPRAGE_reorient.nii.gz"
REORIENTED_TE1="${WORKDIR}/te1_reorient.nii.gz"
REORIENTED_TE2="${WORKDIR}/te2_reorient.nii.gz"

fslreorient2std "${MPRAGE_NII}" "${REORIENTED_MPRAGE}"
fslreorient2std "${TE1_IN}" "${REORIENTED_TE1}"
fslreorient2std "${TE2_IN}" "${REORIENTED_TE2}"

# Segment MPRAGE (T1 anatomical)
echo "Segmenting MPRAGE image (3 tissue classes, T1-type)..."
fast -n 3 -t 1 -o "${WORKDIR}/fast_mprage" "${MPRAGE_NII}"

if [ ! -f "${WORKDIR}/fast_mprage_seg.nii.gz" ]; then
  echo "Error: FAST did not produce a segmentation file for MPRAGE" >&2
  exit 1
fi

echo "Copying MPRAGE FAST outputs to ${OUTDIR}..."
for f in "${WORKDIR}"/fast_mprage*; do
  [ -e "$f" ] && cp "$f" "${OUTDIR}/"
done

# Bias field correction
echo "Bias field correcting TE1..."
fast -B -o "${WORKDIR}/fast_te1" "${REORIENTED_TE1}"
if [ -f "${WORKDIR}/fast_te1_restore.nii.gz" ]; then
  cp "${WORKDIR}/fast_te1_restore.nii.gz" "${TE1_BFC}"
else
  echo "Error: FAST did not produce a restore file for TE1" >&2
  exit 1
fi

echo "Bias field correcting TE2..."
fast -B -o "${WORKDIR}/fast_te2" "${REORIENTED_TE2}"
if [ -f "${WORKDIR}/fast_te2_restore.nii.gz" ]; then
  cp "${WORKDIR}/fast_te2_restore.nii.gz" "${TE2_BFC}"
else
  echo "Error: FAST did not produce a restore file for TE2" >&2
  exit 1
fi

# Registration using bias-corrected images
TE1_TO_MPRAGE="${OUTDIR}/te1_to_MPRAGE.nii.gz"
TE2_TO_MPRAGE="${OUTDIR}/te2_to_MPRAGE.nii.gz"
TE2_TO_MPRAGE_TO_TE1="${OUTDIR}/te2_to_MPRAGE_to_te1.nii.gz"

echo "Registering TE1 (bias corrected) to MPRAGE..."
flirt -in "${TE1_BFC}" -ref "${REORIENTED_MPRAGE}" -out "${TE1_TO_MPRAGE}"

echo "Registering TE2 (bias corrected) to MPRAGE..."
flirt -in "${TE2_BFC}" -ref "${REORIENTED_MPRAGE}" -out "${TE2_TO_MPRAGE}"

echo "Registering TE2_to_MPRAGE to TE1_to_MPRAGE..."
flirt -in "${TE2_TO_MPRAGE}" -ref "${TE1_TO_MPRAGE}" -out "${TE2_TO_MPRAGE_TO_TE1}"

# Normalize TE2 image and compute uT2 map
TE2_NORMALIZED="${OUTDIR}/te2_to_MPRAGE_to_te1_normalized.nii.gz"
UT2_MAP="${OUTDIR}/ut2_map.nii.gz"

fslmaths "${TE2_TO_MPRAGE_TO_TE1}" -mul "${TE2_NORMALIZATION_FACTOR}" "${TE2_NORMALIZED}"
fslmaths "${TE1_TO_MPRAGE}" -sub "${TE2_NORMALIZED}" -div "${TE1_TO_MPRAGE}" "${UT2_MAP}"

echo "Final outputs written to ${OUTDIR}:"
ls -1 "${OUTDIR}"/*.nii.gz
echo "Done."
