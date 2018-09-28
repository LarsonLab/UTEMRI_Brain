%% scripts for dicom generation
clear all
close all
clc

%% load data
% current directory under UTE_Brain/ImageResconstruction/, or
% UTE_Brain/Fitting/


addpath(genpath('../'));


load ../../test_data/brainmask.mat
load ../../test_data/P26112.7-csreconallec_l2_r0p01.mat  % UTE MRI reconstructed images

im_ute = imall(:,:,:,1);

% generate im_biasfield
im_biasfield = ute_brain_estimate_bias_field(im_ute, brainmask);


pfile_name = '../../test_data/P26112.7';
fit_filename = '../../test_data/fit_results.mat';
fit_plus_filename = '../../test_data/fit_results_plus.mat'; 
AIC_thresh = -145;

% generate fits_map
fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield);

% write dicom
write_dicom(fit_maps, pfile_name,im_ute)



