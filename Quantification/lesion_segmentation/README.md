# Lesion Segmentation Pipeline

Pipeline for automated lesion segmentation from MPRAGE and FLAIR MRI scans using **FreeSurfer SAMSEG**. The pipeline supports both NIfTI and DICOM inputs and can process individual subjects or batches of subjects.

## Features

* Accepts NIfTI (`.nii`, `.nii.gz`) or DICOM series directories
* Registers FLAIR to MPRAGE
* Runs FreeSurfer SAMSEG lesion segmentation
* Registers outputs to a user-supplied reference image
* Generates a binary lesion mask in reference-image space
* Supports batch processing of multiple subjects

## Requirements

* FSL
* FreeSurfer (tested with v7.4.1)
* dcm2niix

## Usage

### Single Subject

```bash
bash lesion_segmentation.sh \
    <mprage> \
    <flair> \
    <registration_ref> \
    <output_dir>
```

### Batch Processing

```bash
bash lesion_segmentation_batch.sh \
    <subject_list.txt> \
    lesion_segmentation.sh
```

## Subject List Format

Each subject is defined by four consecutive lines:

```text
<mprage>
<flair>
<registration_ref>
<output_dir>
```

Separate subjects with one or more blank lines. Lines beginning with `#` are ignored.

Example:

```text
# Subject 001
/data/sub001/mprage
/data/sub001/flair
/refs/ref.nii.gz
/results/sub001

# Subject 002
/data/sub002/mprage.nii.gz
/data/sub002/flair.nii.gz
/refs/ref.nii.gz
/results/sub002
```

## Outputs

The pipeline produces:

* `flair_to_mprage.nii.gz`
* `mprage_to_registration_ref.nii.gz`
* `flair_to_mprage_to_registration_ref.nii.gz`
* `seg.nii.gz`
* `seg_to_registration_ref.nii.gz`
* `lesion_mask_to_registration_ref.nii.gz`

Registration transform matrices are also saved for reproducibility.
