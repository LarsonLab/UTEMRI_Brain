#!/bin/tcsh

mkdir ./OUT
rm ./OUT/*

# initial conversion of scanner files from DICOM to Nifti files
dcm2niix -o /OUT /path/to/DICOM/folder

echo "now running BET protocol"
bet /OUT/path/to/input/nifti extracted_data -f 0.3 -g 0 -m
