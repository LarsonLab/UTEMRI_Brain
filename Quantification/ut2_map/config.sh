#!/usr/bin/env bash

# Central configuration file for the UT2 mapping pipeline.
#
# Edit this file once after installation to match your local
# computing environment and software locations.
#
# Variables:
#   FSL_DIR          -> module name/path for the FSL installation
#   APPTAINER_IMAGE  -> path to the MindGlide .sif container image
#   MODEL_BIND       -> bind mount mapping for MindGlide model files

export FSL_DIR="SCS/fsl/fsl_latest"
export APPTAINER_IMAGE="/data/larson9/mindGlide/mindglide.sif"
export MODEL_BIND="/data/larson9/mindGlide/models:/models"