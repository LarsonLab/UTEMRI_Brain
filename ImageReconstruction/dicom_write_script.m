%% scripts for sample generation

%% load data

addpath(genpath('../'));

pfile_name = '../../test_data/P26112.7';

fit_filename = '../../test_data/fit_results.mat';
fit_plus_filename = '../../test_data/fit_results_plus.mat'; 
AIC_thresh = -145;
im_biasfield = [];


fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield);

write_dicom(fit_maps, pfile_name)



