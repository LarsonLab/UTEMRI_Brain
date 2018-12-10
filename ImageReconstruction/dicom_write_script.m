%% scripts for dicom generation
clear all
close all
clc

%% load data
% current directory under UTE_Brain/ImageResconstruction/, or
% UTE_Brain/Fitting/


addpath(genpath('../'));


load /Users/nikhil/Documents/ute_t2/brainmask.mat
load /Users/nikhil/Documents/ute_t2/P26112.7-csreconallec_l2_r0p01.mat  % UTE MRI reconstructed images

im_ute = imall(:,:,:,1);

% generate im_biasfield
im_biasfield = ute_brain_estimate_bias_field(im_ute, brainmask);


pfile_name = '/Users/nikhil/Documents/ute_t2/P26112.7';
fit_filename = '/Users/nikhil/Documents/ute_t2/fit_results.mat';
fit_plus_filename = '/Users/nikhil/Documents/ute_t2/fit_results_plus.mat'; 
AIC_thresh = -145;

% generate fits_map
fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield);

% write dicom
write_ute_brain_dicom(fit_maps, pfile_name,im_ute)



