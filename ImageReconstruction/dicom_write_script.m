%% scripts for sample generation

%% load data

% choose directory
addpath(genpath('../UTE_Brain/'));

pfile_name = './P26112.7';

fit_filename = './fit_results_plus.mat';
fit_plus_filename = './fit_results.mat'; 
AIC_thresh = -145;
im_biasfield = [];


fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield);

write_dicom(fit_maps, pfile_name)



