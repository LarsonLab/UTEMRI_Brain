clear;

%% Select input and output configurations
write_flag = 1; % 1 - save outputs, 0 - don't save output
load_trajectory_flag = 0; % 1 - load existing trajectory, 0 - generate trajectory

twix_path = '/data/larson9/ReVIVE/ReVIVE_S15/3_months/meas_MID02945_FID19495_fid_rosette_2mm_bSSFP_40k.dat'; % path to twix file
trajectory_path = ''; % path to .mat trajectory file containing trajectory_te1 and trajectory_te2

out_path = '/home/abechtel/Documents/UTERosette_bSSFP_sandbox/sandbox/data/ReVIVE_S15_test_subject/3_months/trajectory_ari/test3/'; % path to save output files to

%% Set up Dependencies
mirt_path = '/home/plarson/matlab/reconstruction/mirt/'; % path to mirt reconstruction toolbox

if (~exist('Gmri','file'))
    %set up IRT 
    currentDirectory = pwd;
    cd(mirt_path);
    setup
    cd(currentDirectory)
end

%% Read in raw data
[refmrprot, refmdh, inref] =rdMeas_dene(twix_path);

%% Set scan and recon parameters
recon_parameters.sampling_interval=10; %in us
recon_parameters.fov=0.24; %field of view in m
recon_parameters.matrix_size=240; % number of samples along one spatial dimension
recon_parameters.npoints_skip_te1 = 2; % number of data points to skip in first echo
recon_parameters.npoints_skip_te2 = 0; % number of data points to skip in second echo

% reconstructed matrix size
nx=256;
ny=256;
nz=256;

recon_parameters.ncoils = size(inref,3);
recon_parameters.npetals = size(inref,2);
recon_parameters.nsamples_per_petal = size(inref,1);

% calculate starting data point of each echo (along each petal in rosette)
recon_parameters.start_te1 = recon_parameters.npoints_skip_te1 + 1;
recon_parameters.start_te2 = recon_parameters.nsamples_per_petal/2 + recon_parameters.npoints_skip_te2 + 1;

%% Data modulations
% phase correction for chopping(?)
inref(:,2:2:end,:) = -inref(:,2:2:end,:);

% demodulate k-space data to account for off-resonance
recon_parameters.frequency_offset = 300; %Hz
t=(0:recon_parameters.nsamples_per_petal-1) * recon_parameters.sampling_interval/2;  % us
recon_parameters.f_modulation = exp(1j*2*pi*recon_parameters.frequency_offset*t'/1e6);
inref = inref .* repmat(recon_parameters.f_modulation, [1, recon_parameters.npetals, recon_parameters.ncoils]);

%% Get k-space trajectories

if load_trajectory_flag
    load(trajectory_path, 'trajectory_te1', 'trajectory_te2');
    recon_parameters.trajectory_te1 = trajectory_te1;
    recon_parameters.trajectory_te2 = trajectory_te2;
else
    [recon_parameters.trajectory_te1, recon_parameters.trajectory_te2] = ...
        generate_dual_echo_rosette_trajectory(recon_parameters.matrix_size, recon_parameters.fov);
end

%% Reconstruct first echo
trajectory = recon_parameters.trajectory_te1;
data=squeeze(inref(recon_parameters.start_te1:recon_parameters.start_te1 + recon_parameters.nsamples_per_petal/2 - 1,:,:));
img_recon_te1 = reconstruct_coil_compressed(data, trajectory, nx, ny,nz);

%%  Reconstruct second echo
trajectory = recon_parameters.trajectory_te2;
data=squeeze(inref(recon_parameters.start_te2:recon_parameters.start_te2 + recon_parameters.nsamples_per_petal/2 - 1,:,:));
img_recon_te2 = reconstruct_coil_compressed(data, trajectory, nx, ny,nz);

%% Save outputs
[~, twix_filename, file_ext] = fileparts(twix_path);
recon_parameters.twix_filename = [twix_filename, file_ext];

if write_flag
    save([out_path 'img_recon.mat'],'img_recon_te1', 'img_recon_te2');
    niftiwrite(img_recon_te1, [out_path 'img_recon_te1.nii']);
    niftiwrite(img_recon_te2, [out_path 'img_recon_te2.nii']);
    save([out_path 'recon_parameters.mat'], 'recon_parameters');
end