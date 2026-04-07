# UTE Rosette Dual-Echo MRI Reconstruction

## Overview

This repository implements a MATLAB-based reconstruction pipeline for
**dual-echo UTE rosette MRI data** acquired on Siemens scanners.

The pipeline: 
- Reads raw Siemens **twix data** and applies preprocessing
(phase correction, demodulation, coil compression)
- Generates or loads **3D rosette trajectories**
- Performs **NUFFT-based reconstruction**
- Estimates **coil sensitivity maps**
- Combines coils into final 3D
images
- Outputs both `.mat` and **NIfTI volumes**

The pipeline supports single-subject reconstruction and batch
reconstruction using config files

------------------------------------------------------------------------

## Pipeline Structure

### Main Entry Points

-   `recon_single_scan.m`: Single scan reconstruction script
-   `batch_recon.m`: Batch reconstruction script
-   `UTERosette_bSSFP_dual_echo_reconstruction_Siemens.m`: Core
    reconstruction function

### Configuration

-   `base_config.m`
    Defines default reconstruction parameters, including:
    -   I/O paths
    -   Trajectory parameters
    -   Reconstruction settings (FOV, matrix size, interpolation, etc.)

### Supporting Modules

-   `generate_dual_echo_rosette_trajectory.m`: Trajectory generation
-   `coilwise_nufft_recon.m`: NUFFT reconstruction (per coil)
-   `estimate_csm.m`: Coil sensitivity estimation
-   `combine_coils.m`: Coil combination
### Batch Utilities

-   `generate_config_file_batch.m`
    Generates per-scan config files and a batch list

------------------------------------------------------------------------

## Reconstruction Workflow

1.  Load raw Siemens twix data
2.  Apply preprocessing (phase correction, demodulation)
3.  Load or generate trajectories
4.  Perform coil compression (SVD)
5.  Upsample/interpolate k-space data
6.  Separate TE1 and TE2
7.  Perform NUFFT reconstruction (per coil)
8.  Estimate coil sensitivity maps (from TE1)
9.  Combine coils into final images
10. Save outputs (MAT, NIfTI)

------------------------------------------------------------------------

## Dependencies

- MATLAB  
- MIRT toolbox (NUFFT)  
- Siemens twix reader (`rdMeas_dene`)  

**Important:** This repository contains hardcoded paths to external dependencies. You must update these paths in the following files before running the pipeline:

- `recon_single_scan.m`
- `batch_recon.m`

Specifically, update:

```matlab
rdMeas_dene_path = '/path/to/ImageReconstruction';
cd('/path/to/mirt/');
```

------------------------------------------------------------------------

## Usage

### Single-Subject Reconstruction

Edit and run:

``` matlab
recon_single_scan.m
```

Key steps:

``` matlab
config = base_config();
config.io.twix_path = '/path/to/meas.dat';
config.io.out_path  = '/path/to/output/';
```

Optional:

``` matlab
config.io.trajectory_path = '/path/to/trajectory.mat';
config.io.dcf_path        = '/path/to/dcf.mat';
```

Run:

``` matlab
[img_te1, img_te2] = ...
    UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config);
```

------------------------------------------------------------------------

### Batch Reconstruction

#### Step 1: Generate configs

Edit:

``` matlab
generate_config_file_batch.m
```

Define scans:

``` matlab
scans = {
    '/path/to/meas1.dat', '/output1';
    '/path/to/meas2.dat', '/output2';
};
```

Run script to create:
- `config.mat` files 
- `.txt` file containing list of paths to config files

#### Step 2: Run batch

``` matlab
batch_recon('/path/to/config_list.txt');
```

------------------------------------------------------------------------

## Configuration

Controlled via `config` struct:

### I/O

-   `config.io.twix_path`
-   `config.io.out_path`
-   `config.io.trajectory_path`
-   `config.io.dcf_path`

### Reconstruction

-   FOV, matrix size
-   Sampling interval
-   Frequency offset
-   Upsampling/interpolation
-   Echo skip/shift parameters

### Trajectory

-   Rosette sampling parameters
-   Angular sampling
-   Echo start positions

------------------------------------------------------------------------

## Outputs

    <output_dir>/UTERosette_bSSFP_dual_echo/
        config.mat
        img_recon.mat
        trajectory.mat
        img_recon_te1.nii
        img_recon_te2.nii

------------------------------------------------------------------------

## Notes

-   Coil sensitivity maps are estimated from TE1 and reused for TE2
-   Supports precomputed trajectories and DCF weights
-   Uses MATLAB parallel toolbox (`parpool`)
-   Batch mode runs serially
