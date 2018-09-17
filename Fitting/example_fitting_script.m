clear variables;
writedata_flag = 0;

%%
Ifolder = 1;
foldername{1} = '../SampleData/';  % must have git lfs ...
foldername{1} = '/data/larson/brain_uT2/2018-07-20-3T-DTI-volunteer/';
B0 = 3;

%% load data
%load([foldername{Ifolder} '/ute_32echo_random-csreconallec_l2_r0p01.mat'])
load([foldername{Ifolder} '/P26112.7-csreconallec_l2_r0p01.mat'])

%% parameter estimate
phi_RF = .0562;

%% determine which voxels to fit
load([foldername{Ifolder} '/brainmask.mat'])


fit_mask = repmat(brainmask, [1,1,1,size(imall,4)]);
% on -resonance fit
savefname = [foldername{Ifolder} '/fit_results'];
ute_brain_fitting(TE, imall .* fit_mask, B0, phi_RF, savefname);

% off-resonance ("plus" frequency) fit
savefname = [foldername{Ifolder} '/fit_results_plus'];
ute_brain_fitting(TE, imallplus .* fit_mask, B0, phi_RF, savefname);





