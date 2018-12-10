#!/bin/tcsh

mkdir ./OUT
rm ./OUT/*

# initial conversion of scanner files from DICOM to Nifti files
dcm2niix -o /Users/nikhil/Documents/ute_t2/E6674

echo "now running BET protocol"
bet /Users/nikhil/Documents/ute_t2/E6674/1/E6674S1I1.DCM extracted_data -f 0.3 -g 0 -m
