#!/bin/bash
set -euo pipefail

# Script to generate uT2 map from MPRAGE and dual-echo UTE images
# Usage: ./generate_ut2_map.sh -m <MPRAGE_dicom_or_nifti> -te1 <te1.nii.gz> -te2 <te2.nii.gz> -o <output_directory>

#!/bin/bash
set -euo pipefail

# Load FSL module (currently skipping this for batch processing purposes)
# FSL_DIR="SCS/fsl/fsl_latest"
# module load "$FSL_DIR"

# Initialize variables
MPRAGE_INPUT=""
TE1_IN=""
TE2_IN=""
OUTDIR=""

# Parse command line flags
while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--mprage) MPRAGE_INPUT="$2"; shift 2 ;;
    -te1|--te1) TE1_IN="$2"; shift 2 ;;
    -te2|--te2) TE2_IN="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -m <MPRAGE_dicom_or_nifti> -te1 <te1.nii.gz> -te2 <te2.nii.gz> -o <output_directory>"
      exit 0 ;;
    *) echo "Unknown option: $1"; echo "Use -h for help."; exit 1 ;;
  esac
done

# Validate inputs
if [[ -z "$MPRAGE_INPUT" || -z "$TE1_IN" || -z "$TE2_IN" || -z "$OUTDIR" ]]; then
  echo "ERROR: Missing required arguments."
  echo "Usage: $0 -m <MPRAGE_dicom_or_nifti> -te1 <te1.nii.gz> -te2 <te2.nii.gz> -o <output_directory>"
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
TE1_BFC="${OUTDIR}/te1_bfc.nii.gz"
TE2_BFC="${OUTDIR}/te2_bfc.nii.gz"
TE1_BIASMAP="${OUTDIR}/te1_bfcmap.nii.gz"
TE2_BIASMAP="${OUTDIR}/te2_bfcmap.nii.gz"

# Detect MPRAGE input type
if [ -d "${MPRAGE_INPUT}" ]; then
  echo "Detected MPRAGE input as DICOM directory."
  dcm2niix -z y -f MPRAGE -o "${WORKDIR}" "${MPRAGE_INPUT}"
  MPRAGE_NII="$(ls -1 "${WORKDIR}"/*.nii.gz 2>/dev/null | head -n1 || true)"
  if [ -z "${MPRAGE_NII}" ]; then
    echo "ERROR: dcm2niix did not produce a .nii.gz file in ${WORKDIR}."
    exit 1
  fi
  cp "${MPRAGE_NII}" "${MPRAGE_BASE}.nii.gz"
  MPRAGE_NII="${MPRAGE_BASE}.nii.gz"
elif [[ "${MPRAGE_INPUT}" == *.nii ]] || [[ "${MPRAGE_INPUT}" == *.nii.gz ]]; then
  echo "Detected MPRAGE input as NIfTI file."
  # Extract original extension
  EXT="${MPRAGE_INPUT##*.}"        # gets 'nii' or 'gz'
  if [[ "${EXT}" == "gz" ]]; then
    cp "${MPRAGE_INPUT}" "${MPRAGE_BASE}.nii.gz"
    MPRAGE_NII="${MPRAGE_BASE}.nii.gz"
  else
    cp "${MPRAGE_INPUT}" "${WORKDIR}/MPRAGE.nii"
    fslchfiletype NIFTI_GZ "${WORKDIR}/MPRAGE.nii" "${MPRAGE_BASE}"
    MPRAGE_NII="${MPRAGE_BASE}.nii.gz"
  fi
else
  echo "ERROR: -m argument must be either a DICOM folder or a .nii/.nii.gz file."
  exit 1
fi

# Reorient to standard
echo "Reorienting images..."
REORIENTED_MPRAGE="${WORKDIR}/MPRAGE_reorient.nii.gz"
REORIENTED_TE1="${WORKDIR}/te1_reorient.nii.gz"
REORIENTED_TE2="${WORKDIR}/te2_reorient.nii.gz"
fslreorient2std "${MPRAGE_NII}" "${REORIENTED_MPRAGE}"
fslreorient2std "${TE1_IN}" "${REORIENTED_TE1}"
fslreorient2std "${TE2_IN}" "${REORIENTED_TE2}"

# Segment MPRAGE
echo "Segmenting MPRAGE..."
fast -n 3 -t 1 -o "${WORKDIR}/fast_mprage" "${MPRAGE_NII}"
if [ ! -f "${WORKDIR}/fast_mprage_seg.nii.gz" ]; then
  echo "ERROR: FAST segmentation failed for MPRAGE."
  exit 1
fi
cp "${WORKDIR}"/fast_mprage* "${OUTDIR}/"

# Bias field correction for TE1
echo "Bias field correcting TE1..."
fast -B -b -o "${WORKDIR}/fast_te1" "${REORIENTED_TE1}"
cp "${WORKDIR}/fast_te1_restore.nii.gz" "${TE1_BFC}"
cp "${WORKDIR}/fast_te1_bias.nii.gz" "${TE1_BIASMAP}"

# Bias field correction for TE2
echo "Bias field correcting TE2..."
fast -B -b -o "${WORKDIR}/fast_te2" "${REORIENTED_TE2}"
cp "${WORKDIR}/fast_te2_restore.nii.gz" "${TE2_BFC}"
cp "${WORKDIR}/fast_te2_bias.nii.gz" "${TE2_BIASMAP}"

# Registration
TE1_TO_MPRAGE="${OUTDIR}/te1_to_MPRAGE.nii.gz"
TE2_TO_MPRAGE="${OUTDIR}/te2_to_MPRAGE.nii.gz"
TE2_TO_MPRAGE_TO_TE1="${OUTDIR}/te2_to_MPRAGE_to_te1.nii.gz"

echo "Registering TE1 to MPRAGE..."
flirt -in "${TE1_BFC}" -ref "${REORIENTED_MPRAGE}" -out "${TE1_TO_MPRAGE}"

echo "Registering TE2 to MPRAGE..."
flirt -in "${TE2_BFC}" -ref "${REORIENTED_MPRAGE}" -out "${TE2_TO_MPRAGE}"

echo "Registering TE2_to_MPRAGE to TE1_to_MPRAGE..."
flirt -in "${TE2_TO_MPRAGE}" -ref "${TE1_TO_MPRAGE}" -out "${TE2_TO_MPRAGE_TO_TE1}"

# Compute TE2 normalization factor
MPRAGE_SEG="${OUTDIR}/fast_mprage_seg.nii.gz"
GM_MASK="${WORKDIR}/mprage_gm_mask.nii.gz"
fslmaths "${MPRAGE_SEG}" -thr 2 -uthr 2 -bin "${GM_MASK}"

MEAN_TE1=$(fslstats "${TE1_TO_MPRAGE}" -k "${GM_MASK}" -M)
MEAN_TE2=$(fslstats "${TE2_TO_MPRAGE_TO_TE1}" -k "${GM_MASK}" -M)
if (( $(echo "${MEAN_TE2} == 0" | bc -l) )); then
  echo "Error: Mean TE2 signal is zero."
  exit 1
fi
TE2_NORMALIZATION_FACTOR=$(echo "${MEAN_TE1} / ${MEAN_TE2}" | bc -l)
echo "TE2 normalization factor: ${TE2_NORMALIZATION_FACTOR}"

# Normalize and compute uT2 map
TE2_NORMALIZED="${OUTDIR}/te2_to_MPRAGE_to_te1_normalized.nii.gz"
UT2_MAP="${OUTDIR}/ut2_map.nii.gz"
fslmaths "${TE2_TO_MPRAGE_TO_TE1}" -mul "${TE2_NORMALIZATION_FACTOR}" "${TE2_NORMALIZED}"
fslmaths "${TE1_TO_MPRAGE}" -sub "${TE2_NORMALIZED}" -div "${TE1_TO_MPRAGE}" "${UT2_MAP}"

# Final report
echo
echo "Final outputs written to ${OUTDIR}:"
ls -1 "${OUTDIR}"/*.nii.gz
echo "Done."

