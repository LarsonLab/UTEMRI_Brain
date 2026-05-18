#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_mindglide.sh /path/to/input.nii.gz /path/to/output_seg.nii.gz

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.sh"

APPTAINER_IMAGE="/home/$USER/mindGlide/mindglide.sif"
MODEL_BIND="/home/$USER/mindGlide/models:/models"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/input_image.nii.gz /path/to/output_segmentation.nii.gz"
    exit 1
fi

INPUT_ARG="$1"
OUTPUT_ARG="$2"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ ! -f "$INPUT_ARG" ]]; then
    echo "ERROR: input file not found: $INPUT_ARG"
    exit 1
fi

# Resolve absolute paths
INPUT_IMAGE="$(cd "$(dirname "$INPUT_ARG")" && pwd)/$(basename "$INPUT_ARG")"
OUTPUT_DIR="$(dirname "$OUTPUT_ARG")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_SEG="$(cd "$OUTPUT_DIR" && pwd)/$(basename "$OUTPUT_ARG")"

INPUT_DIR="$(dirname "$INPUT_IMAGE")"
OUTPUT_DIR_ABS="$(dirname "$OUTPUT_SEG")"

apptainer exec --nv \
    -B "${MODEL_BIND}" \
    -B "${INPUT_DIR}:${INPUT_DIR}" \
    -B "${OUTPUT_DIR_ABS}:${OUTPUT_DIR_ABS}" \
    "${APPTAINER_IMAGE}" \
    bash -lc "export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; \
      mindglide -i '${INPUT_IMAGE}' -o '${OUTPUT_SEG}' --sw_batch_size 1"