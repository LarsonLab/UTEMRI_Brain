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
TE1_TO_MPRAGE="${OUTDIR}/te1_to_MPRAGE"
TE2_TO_MPRAGE="${OUTDIR}/te2_to_MPRAGE"
TE2_TO_MPRAGE_TO_TE1="${OUTDIR}/te2_to_MPRAGE_to_te1"

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

MPRAGE_NII="${REORIENTED_MPRAGE}"
TE1_IN="${REORIENTED_TE1}"
TE2_IN="${REORIENTED_TE2}"

# --- Registration ---
echo "Registering te1 to MPRAGE..."
flirt -in "${TE1_IN}" -ref "${MPRAGE_NII}" -out "${TE1_TO_MPRAGE}.nii.gz"

echo "Registering te2 to MPRAGE..."
flirt -in "${TE2_IN}" -ref "${MPRAGE_NII}" -out "${TE2_TO_MPRAGE}.nii.gz"

echo "Registering te2_to_MPRAGE to te1_to_MPRAGE..."
flirt -in "${TE2_TO_MPRAGE}.nii.gz" -ref "${TE1_TO_MPRAGE}.nii.gz" -out "${TE2_TO_MPRAGE_TO_TE1}.nii.gz"

# Bias field correction using FAST
echo "Bias field correcting te1_to_MPRAGE..."
fast -B -o "${WORKDIR}/fast_te1" "${TE1_TO_MPRAGE}.nii.gz"
if [ -f "${WORKDIR}/fast_te1_restore.nii.gz" ]; then
  TE1_RESTORE="${WORKDIR}/fast_te1_restore.nii.gz"
else
  echo "Error: FAST did not produce a restore file for TE1" >&2
  exit 1
fi
cp "${TE1_RESTORE}" "${OUTDIR}/te1_to_MPRAGE_restore.nii.gz"

echo "Bias field correcting te2_to_MPRAGE_to_te1..."
fast -B -o "${WORKDIR}/fast_te2" "${TE2_TO_MPRAGE_TO_TE1}.nii.gz"
if [ -f "${WORKDIR}/fast_te2_restore.nii.gz" ]; then
  TE2_RESTORE="${WORKDIR}/fast_te2_restore.nii.gz"
else
  echo "Error: FAST did not produce a restore file for TE2" >&2
  exit 1
fi
cp "${TE2_RESTORE}" "${OUTDIR}/te2_to_MPRAGE_to_te1_restore.nii.gz"

# Normalize TE2 image and compute uT2 map
TE2_NORMALIZED="${OUTDIR}/te2_to_MPRAGE_to_te1_restore_normalized.nii.gz"
UT2_MAP="${OUTDIR}/ut2_map.nii.gz"

fslmaths "${OUTDIR}/te2_to_MPRAGE_to_te1_restore.nii.gz" -mul "${TE2_NORMALIZATION_FACTOR}" "${TE2_NORMALIZED}"
fslmaths "${OUTDIR}/te1_to_MPRAGE_restore.nii.gz" -sub "${TE2_NORMALIZED}" "${UT2_MAP}"

echo "Final outputs written to ${OUTDIR}:"
ls -1 "${OUTDIR}"/*.nii.gz
echo "Done."
