% Batch reconstruction script for UTE Rosette dual-echo Siemens data
%
% This script loops through a list of raw twix files, creates a config
% struct for each scan, and calls:
%   UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config)
%
% You can run this serially or with a parfor loop if you want parallelism.

% List of raw twix files and output directories
scans = {
    '/data/larson9/ReVIVE/ReVIVE_S17/3_months/meas_MID00336_FID27805_fid_rosette_2mm_bSSFP_40k.dat', '/data/larson9/ReVIVE/ReVIVE_S17/3_months/';
};

% Modify reconstruction parameters if needed
config = base_config();
% config.recon.npoints_skip_te1 = 2;
% config.recon.npoints_skip_te2 = 0;

% Open parallel pool (adjust worker count to SLURM allocation)
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('threads'); % or parpool(N) if you know worker count
end

% Setup dependencies
if (~exist('Gmri','file'))
    %set up mirt recon toolbox
    currentDirectory = pwd;
    cd('/home/plarson/matlab/reconstruction/mirt/'); setup
    cd(currentDirectory)
end

% Run reconstruction in parallel
nscans = size(scans,1);
input = scans(:,1);
outdir = scans(:,2);
parfor i = 1:nscans
    subj_config = config;
    subj_config.io.twix_path = input{i};
    subj_config.io.out_path  = outdir{i};
    fprintf('Reconstructing scan %d/%d: %s\n', i, nscans, subj_config.io.twix_path);

    [~,~,~,~] = UTERosette_bSSFP_dual_echo_reconstruction_Siemens(subj_config);
end

fprintf('Batch reconstruction complete.\n');
