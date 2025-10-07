% Script to reconstruct images for single scan dual echo rosette UTE

clear;

% Initialize default configuration settings
config = base_config();

% Set path to raw Siemens twix file
config.io.twix_path = '/path/to/raw/twix/data/meas_MID123.dat';

% Set output directory
config.io.out_path = '/path/to/output/';
if ~exist(config.io.out_path, 'dir')
    mkdir(config.io.out_path);
end

% Modify recon parameters (optional) - see base_config.m for options
config.recon.npoints_skip_te1 = 2; % number of data points to skip in first echo
config.recon.npoints_skip_te2 = 0; % number of data points to skip in second echo

% Run reconstruction
fprintf('Starting reconstruction for subject: %s\n', config.io.twix_path);
[img_recon_te1, img_recon_te2, trajectory_te1, trajectory_te2] = ...
    UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config);
fprintf('Reconstruction complete. Results saved to: %s\n', config.io.out_path);