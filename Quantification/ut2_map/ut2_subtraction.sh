#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.sh"
module load "$FSL_DIR"

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

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

log "Starting UT2 subtraction"
log "TE1: $TE1_INPUT"
log "TE2: $TE2_INPUT"
log "Rescaling mask: $RESCALING_MASK"

# Measure the mean signal in the rescaling mask for both images.
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

# Scaling factor is mean(TE1) / mean(TE2) within the mask.
TE2_SCALING_FACTOR=$(awk -v a="$TE1_MEAN" -v b="$TE2_MEAN" 'BEGIN { printf "%.10f", a / b }')
if [[ -z "$TE2_SCALING_FACTOR" ]]; then
  echo "Error: failed to compute te2 scaling factor" >&2
  exit 1
fi

log "Mean TE1 in mask: $TE1_MEAN"
log "Mean TE2 in mask: $TE2_MEAN"
log "te2_scaling_factor=$TE2_SCALING_FACTOR"

echo "$TE2_SCALING_FACTOR" > "$OUTDIR/te2_scaling_factor.txt"

# Rescale TE2, then subtract it from TE1 to create the UT2 map.
TE2_RESCALED="$OUTDIR/te2_rescaled.nii.gz"
log "Writing rescaled TE2: $TE2_RESCALED"
fslmaths "$TE2_INPUT" -mul "$TE2_SCALING_FACTOR" "$TE2_RESCALED"

UT2_MAP="$OUTDIR/ut2_map.nii.gz"
log "Writing UT2 map: $UT2_MAP"
fslmaths "$TE1_INPUT" -sub "$TE2_RESCALED" "$UT2_MAP"

log "UT2 subtraction step completed"
log "Outputs written to: $OUTDIR"
