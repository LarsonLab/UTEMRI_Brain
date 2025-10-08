#! /bin/bash
set -euo pipefail

# Constants
TE2_NORMALIZATION_FACTOR=1.38 # scaling factor to normalize TE2 image

# Load FSL module
FSL_DIR="SCS/fsl/fsl_latest" # set to your FSL module path
module load $FSL_DIR

# Initialize empty variables
MPRAGE_DICOM_FOLDER=""
TE1_IN=""
TE2_IN=""
OUTDIR=""

# Parse command line flags
while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--mprage)
      MPRAGE_DICOM_FOLDER="$2"
      shift 2
      ;;
    -te1|--te1)
      TE1_IN="$2"
      shift 2
      ;;
    -te2|--te2)
      TE2_IN="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 -m <MPRAGE_dicom_folder> -te1 <te1.nii> -te2 <te2.nii> -o <output_directory>"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use -h for help."
      exit 1
      ;;
  esac
done

# Validate inputs
if [[ -z "$MPRAGE_DICOM_FOLDER" || -z "$TE1_IN" || -z "$TE2_IN" || -z "$OUTDIR" ]]; then
  echo "ERROR: Missing required arguments."
  echo "Usage: $0 -m <MPRAGE_dicom_folder> -te1 <te1.nii> -te2 <te2.nii> -o <output_directory>"
  exit 1
fi

# Set up directories
mkdir -p "${OUTDIR}" # create output directory if it doesn't exist
WORKDIR="$(mktemp -d)" # create temporary directory for intermediate outputs
trap 'rm -rf "${WORKDIR}"' EXIT # delete temporary directory after script finishes

echo "Working directory: ${WORKDIR}"
echo "Output directory: ${OUTDIR}"

# Set file naming conventions
MPRAGE_BASE="${OUTDIR}/MPRAGE"
TE1_TO_MPRAGE="${OUTDIR}/te1_to_MPRAGE"
TE2_TO_MPRAGE="${OUTDIR}/te2_to_MPRAGE"
TE2_TO_MPRAGE_TO_TE1="${OUTDIR}/te2_to_MPRAGE_to_te1"

# Convert MPRAGE from DICOMs to NIfTI
echo "Converting MPRAGE DICOMs to NIfTI..."
dcm2niix -z n -f MPRAGE -o "${WORKDIR}" "${MPRAGE_DICOM_FOLDER}"
MPRAGE_NII="$(ls -1 "${WORKDIR}"/*.nii 2>/dev/null | head -n1 || true)"
if [ -z "${MPRAGE_NII}" ]; then
  echo "ERROR: dcm2niix did not produce a .nii in ${WORKDIR}."
  exit 1
fi
cp "${MPRAGE_NII}" "${MPRAGE_BASE}.nii"
MPRAGE_NII="${MPRAGE_BASE}.nii"

# Reorient all inputs to standard orientation
echo "Reorienting all input images to standard orientation..."
REORIENTED_MPRAGE="${WORKDIR}/MPRAGE_reorient.nii"
REORIENTED_TE1="${WORKDIR}/te1_reorient.nii"
REORIENTED_TE2="${WORKDIR}/te2_reorient.nii"

fslreorient2std "${MPRAGE_NII}" "${REORIENTED_MPRAGE}"
fslreorient2std "${TE1_IN}" "${REORIENTED_TE1}"
fslreorient2std "${TE2_IN}" "${REORIENTED_TE2}"

# Use reoriented versions for all downstream steps
MPRAGE_NII="${REORIENTED_MPRAGE}"
TE1_IN="${REORIENTED_TE1}"
TE2_IN="${REORIENTED_TE2}"

# Registration steps
echo "Registering te1 to MPRAGE..."
flirt -in "${TE1_IN}" -ref "${MPRAGE_NII}" -out "${TE1_TO_MPRAGE}.nii"

echo "Registering te2 to MPRAGE..."
flirt -in "${TE2_IN}" -ref "${MPRAGE_NII}" -out "${TE2_TO_MPRAGE}.nii"

echo "Registering te2_to_MPRAGE to te1_to_MPRAGE..."
flirt -in "${TE2_TO_MPRAGE}.nii" -ref "${TE1_TO_MPRAGE}.nii" -out "${TE2_TO_MPRAGE_TO_TE1}.nii"

# --- Bias field correction using FSL FAST ---
echo "Bias field correcting te1_to_MPRAGE..."
fast -B -o "${WORKDIR}/fast_te1" "${TE1_TO_MPRAGE}.nii"
if [ -f "${WORKDIR}/fast_te1_restore.nii.gz" ]; then
    TE1_RESTORE_CAND="${WORKDIR}/fast_te1_restore.nii.gz"
elif [ -f "${WORKDIR}/fast_te1_restore.nii" ]; then
    TE1_RESTORE_CAND="${WORKDIR}/fast_te1_restore.nii"
else
    echo "Error: FAST did not produce a restore file for TE1" >&2
    exit 1
fi
fslchfiletype NIFTI "${TE1_RESTORE_CAND}" "${OUTDIR}/te1_to_MPRAGE_restore"
TE1_RESTORE="${OUTDIR}/te1_to_MPRAGE_restore.nii"

echo "Bias field correcting te2_to_MPRAGE_to_te1..."
fast -B -o "${WORKDIR}/fast_te2" "${TE2_TO_MPRAGE_TO_TE1}.nii"
if [ -f "${WORKDIR}/fast_te2_restore.nii.gz" ]; then
    TE2_RESTORE_CAND="${WORKDIR}/fast_te2_restore.nii.gz"
elif [ -f "${WORKDIR}/fast_te2_restore.nii" ]; then
    TE2_RESTORE_CAND="${WORKDIR}/fast_te2_restore.nii"
else
    echo "Error: FAST did not produce a restore file for TE2" >&2
    exit 1
fi
fslchfiletype NIFTI "${TE2_RESTORE_CAND}" "${OUTDIR}/te2_to_MPRAGE_to_te1_restore"
TE2_RESTORE="${OUTDIR}/te2_to_MPRAGE_to_te1_restore.nii"

# --- Normalize TE2 image and compute uT2 map ---
TE2_NORMALIZED="${OUTDIR}/te2_to_MPRAGE_to_te1_restore_normalized.nii"
UT2_MAP="${OUTDIR}/ut2_map.nii"

fslmaths "${TE2_RESTORE}" -mul "${TE2_NORMALIZATION_FACTOR}" "${TE2_NORMALIZED}"
fslmaths "${TE1_RESTORE}" -sub "${TE2_NORMALIZED}" "${UT2_MAP}"

# --- Clean up and enforce consistent format ---
for f in "${OUTDIR}"/*.nii.gz; do
  [ -e "$f" ] || continue
  fslchfiletype NIFTI "$f" "${f%.gz}"
  rm -f "$f"
done

echo "Final outputs written to ${OUTDIR}:"
ls -1 "${OUTDIR}"
echo "Done."