% Script to create a batch of config.mat files for UTE Rosette dual-echo 
% reconstruction and save .txt file of config file paths
%
% Each entry in the "scans" cell array should contain:
%   { twix_file_path , output_directory }

clear;

% List of twix files and corresponding output directories
scans = {
    '/path/to/raw/twix/data/meas_MID123.dat', ...
    '/path/to/output/directory';
};


% Output text file listing all config.mat paths
config_list_txt = '/home/abechtel/Documents/UTERosette_bSSFP/data/config_file_list_test1.txt';


% Optional: shared trajectory/dcf paths
trajectory_path = '/home/abechtel/Documents/UTERosette_bSSFP/data/ari_trajectory_and_dcf/no_truncation/trajectory.mat';   % e.g. '/data/trajectories/rosette_traj.mat'
dcf_path        = '/home/abechtel/Documents/UTERosette_bSSFP/data/ari_trajectory_and_dcf/no_truncation/dcf_weights.mat';   % e.g. '/data/trajectories/rosette_dcf.mat'

% Create config files and list
nscans = size(scans, 1);

% Open text file for writing
fid = fopen(config_list_txt, 'w');
if fid == -1
    error('Could not open config list file for writing: %s', config_list_txt);
end

for i = 1:nscans
    twix_path = scans{i,1};
    out_path  = scans{i,2};

    fprintf('Creating config %d/%d\n', i, nscans);
    fprintf('  Twix: %s\n', twix_path);
    fprintf('  Out : %s\n', out_path);

    % Create output directory if it does not exist
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end

    % Initialize base config
    config = base_config();

    % Populate IO fields
    config.io.twix_path = twix_path;
    config.io.out_path  = out_path;
    config.io.trajectory_path = trajectory_path;
    config.io.dcf_path        = dcf_path;

    % Modify recon parameters (optional) - see base_config.m for options
    config.recon.skip_end_te2 = 2;
    config.recon.shift_te1 = 2*config.recon.upsample_factor + 2;
    config.recon.shift_te2 = 2;

    % Save config
    config_file = fullfile(out_path, 'config.mat');
    save(config_file, 'config');

    % Write config path to text file
    fprintf(fid, '%s\n', config_file);

    fprintf('  Saved: %s\n\n', config_file);
end

fclose(fid);

fprintf('All config files created successfully.\n');
fprintf('Config list written to: %s\n', config_list_txt);
