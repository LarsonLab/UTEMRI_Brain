function batch_recon(config_list_txt)
%BATCH_RECON
% Runs serial batch reconstruction of UTE Rosette dual-echo Siemens data
% using pre-populated config structs stored in .mat files.
%
% Each line of the input text file should contain the full path to a .mat
% file that includes a struct named "config". The config struct must already
% be fully populated with all fields required by:
%   UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config)
%
% INPUTS:
%   config_list_txt : Path to a .txt file where each line is the full path to a
%                   .mat file containing a populated "config" struct.
%
% OUTPUTS:
%   None
%
% NOTES:
%   - Reconstructions are run serially (parallel processing not supported).
%   - Required reconstruction and MIRT dependencies are added to the MATLAB
%     path once at runtime.
%   - Missing or invalid config files are skipped with a warning.

    % Input checks
    if nargin ~= 1
        error('Exactly one input required: path to config list .txt file');
    end

    if ~exist(config_list_txt, 'file')
        error('Config list file not found: %s', config_list_txt);
    end

    % Read list of config .mat files
    fid = fopen(config_list_txt, 'r');
    if fid == -1
        error('Could not open file: %s', config_list_txt);
    end

    config_list = {};
    line = fgetl(fid);
    while ischar(line)
        line = strtrim(line);
        if ~isempty(line)
            config_list{end+1,1} = line;
        end
        line = fgetl(fid);
    end
    fclose(fid);

    nscans = numel(config_list);
    if nscans == 0
        error('No config files found in %s', config_list_txt);
    end

    fprintf('Found %d reconstruction jobs.\n', nscans);

    % Setup dependencies
    rdMeas_dene_path = '/home/abechtel/Documents/UTEMRI_Brain/ImageReconstruction';
    addpath(genpath(rdMeas_dene_path));

    if (~exist('Gmri', 'file'))
        currentDirectory = pwd;
        cd('/home/plarson/matlab/reconstruction/mirt/');
        setup
        cd(currentDirectory)
    end

    % Run reconstructions in serial
    for i = 1:nscans
        subject_config_file = config_list{i};

        if ~exist(subject_config_file, 'file')
            warning('Config file not found, skipping (%d/%d): %s', ...
                    i, nscans, subject_config_file);
            continue;
        end

        fprintf('\n[%d/%d] Loading config: %s\n', i, nscans, subject_config_file);
        S = load(subject_config_file);

        if ~isfield(S, 'config')
            warning('No "config" struct found in %s, skipping.', subject_config_file);
            continue;
        end

        config = S.config;

        % Optional sanity check
        if ~isfield(config, 'io') || ...
           ~isfield(config.io, 'twix_path') || ...
           ~isfield(config.io, 'out_path')
            warning('Config in %s appears incomplete, skipping.', subject_config_file);
            continue;
        end

        fprintf('Reconstructing:\n  Input: %s\n  Output: %s\n', ...
                config.io.twix_path, config.io.out_path);

        try
            [~,~,~,~] = UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config);
        catch ME
            warning('Reconstruction failed for %s\n%s', ...
                    subject_config_file, ME.message);
        end
    end

    fprintf('\nBatch reconstruction complete.\n');
end
