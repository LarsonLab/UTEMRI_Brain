#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") <te1_to_registration_ref.nii.gz> <te2_to_registration_ref_to_te1.nii.gz> <rescaling_mask.nii.gz> <output_dir>
USAGE
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Error: required command not found: $1" >&2; exit 1; }
}

if [[ $# -ne 4 ]]; then
  usage >&2
  exit 1
fi

FSL_DIR="SCS/fsl/fsl_latest"
module load "$FSL_DIR"

TE1_INPUT=$1
TE2_INPUT=$2
RESCALING_MASK=$3
OUTDIR=$4

for f in "$TE1_INPUT" "$TE2_INPUT" "$RESCALING_MASK"; do
  [[ -e "$f" ]] || { echo "Error: input not found: $f" >&2; exit 1; }
done

require_cmd fslstats
require_cmd fslmaths
require_cmd awk

mkdir -p "$OUTDIR"

TE1_MEAN=$(fslstats "$TE1_INPUT" -k "$RESCALING_MASK" -M)
TE2_MEAN=$(fslstats "$TE2_INPUT" -k "$RESCALING_MASK" -M)

if [[ -z "$TE1_MEAN" || -z "$TE2_MEAN" ]]; then
  echo "Error: failed to compute means inside rescaling mask" >&2
  exit 1
fi

if awk -v b="$TE2_MEAN" 'BEGIN { exit !(b == 0) }'; then
  echo "Error: te2 mean inside rescaling mask is zero; cannot compute scaling factor" >&2
  exit 1
fi

TE2_SCALING_FACTOR=$(awk -v a="$TE1_MEAN" -v b="$TE2_MEAN" 'BEGIN { printf "%.10f", a / b }')
if [[ -z "$TE2_SCALING_FACTOR" ]]; then
  echo "Error: failed to compute te2 scaling factor" >&2
  exit 1
fi

echo "te2_scaling_factor=$TE2_SCALING_FACTOR"
echo "$TE2_SCALING_FACTOR" > "$OUTDIR/te2_scaling_factor.txt"

TE2_RESCALED="$OUTDIR/te2_rescaled.nii.gz"
fslmaths "$TE2_INPUT" -mul "$TE2_SCALING_FACTOR" "$TE2_RESCALED"

UT2_MAP="$OUTDIR/ut2_map.nii.gz"
fslmaths "$TE1_INPUT" -sub "$TE2_RESCALED" "$UT2_MAP"

echo "UT2 subtraction step completed. Outputs written to: $OUTDIR"
