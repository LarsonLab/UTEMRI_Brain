function [img_recon_te1, img_recon_te2, ...
          trajectory_te1, trajectory_te2] = ...
          UTERosette_bSSFP_dual_echo_reconstruction_Siemens(config)
% UTEROSETTE_BSSFP_DUAL_ECHO_RECONSTRUCTION_SIEMENS
%
% Reconstructs dual-echo UTE rosette MRI data acquired on a Siemens scanner.
% This function reads in raw Siemens "twix" data, applies preprocessing
% steps (demodulation, phase correction), generates rosette k-space
% trajectories, and reconstructs 3D images for both echo times. The function also saves
% reconstructed images and trajectories.
%
% Inputs:
%   config                - Struct containing input/output and recon configuration:
%       .io.twix_path     : Path to raw Siemens twix data file
%       .io.out_path      : Directory for saving output files
%       .io.trajectory_path : (optional) Path to .mat file containing trajectories
%                             (fields: trajectory_te1, trajectory_te2)
%       .io.dcf_path      : (optional) Path to .mat file containing dcf weights 
%                           (fields: dcf_weights_te1, dcf_weights_te2)
%       .recon.output_size : Reconstructed matrix size
%                           (fields: nx, ny, nz)
%       .recon.matrix_size : Matrix size used for trajectory generation
%       .recon.fov        : Field of view (mm)
%       .recon.npoints_skip_te1 : # of samples to skip at start of TE1
%       .recon.npoints_skip_te2 : # of samples to skip at start of TE2
%       .recon.sampling_interval: Sampling interval (µs)
%       .recon.frequency_offset : Frequency offset used to demodulate
%                                 kspace data (Hz)
%
% Outputs:
%   img_recon_te1    - 3D complex reconstructed image at echo 1
%   img_recon_te2    - 3D complex reconstructed image at echo 2
%   trajectory_te1   - K-space trajectory for first echo
%   trajectory_te2   - K-space trajectory for second echo
%
% Saves the following to config.io.out_path:
%   - 'config.mat' containing configuration struct
%   - 'img_recon.mat' containing reconstructed images for both echoes
%   - 'trajectory.mat' containing k-space trajectories for both echos
%   - 'img_recon_te1.nii' and 'img_recon_te2.nii' NIfTI volumes
%
% Example:
%   config = base_config();
%   config.io.twix_path = 'subject01/meas_MID123.dat';
%   config.io.out_path  = 'subject01/output/';
%   [img1, img2, traj1, traj2] = ...
%       UTERosette_bSSFP_dual_echo_reconstruction_Siemens(cfg);

    % -- Detect number of workers to use in parallel processing (e.g., from SLURM) --
    nCPUs = str2double(getenv('SLURM_CPUS_ON_NODE'));
    if isnan(nCPUs) || nCPUs == 0
        nCPUs = feature('numcores');
    end

    % Start the pool
    if isempty(gcp('nocreate'))
        parpool('local', nCPUs);
    end

    tic;

    % -- Read in raw data --
    [~, ~, raw_data] =rdMeas_dene(config.io.twix_path);
    
    % -- Setup some recon parameters --
    % Reconstructed matrix size
    nx=config.recon.output_size.nx;
    ny=config.recon.output_size.ny;
    nz=config.recon.output_size.nz;
    
    ncoils = size(raw_data,3);
    npetals = size(raw_data,2);
    nsamples_per_petal = size(raw_data,1);
    
    % calculate starting data point of each echo (along each petal in rosette)
    start_te1 = 1;
    start_te2 = nsamples_per_petal/2 + 1;
    
    % -- Data modulations --
    % Phase correction for chopping(?)
    raw_data(:,2:2:end,:) = -raw_data(:,2:2:end,:);
    
    % Demodulate k-space data to account for off-resonance
    t=(0:nsamples_per_petal-1) * config.recon.sampling_interval/2;  % us
    f_modulation = exp(1j*2*pi*config.recon.frequency_offset*t'/1e6);
    raw_data = raw_data .* repmat(f_modulation, [1, npetals, ncoils]);
    
    % -- Load/generate k-space trajectories --

    % Check if valid trajectory file path provided in config file
    if isfile(config.io.trajectory_path)
        % Load trajectories
        disp('Trajectory file detected.');
        disp('Loading input trajectories...');
        load(config.io.trajectory_path, 'trajectory_te1', 'trajectory_te2');
    else
        % Generate trajectories
        disp('No valid trajectory file detected.');
        disp('Generating trajectories...');
        [trajectory_te1, trajectory_te2] = ...
            generate_dual_echo_rosette_trajectory(config.recon.traj);
    end

    % -- Coil compression via SVD --
    disp('Performing coil compression...');
    data_compressed = squeeze(raw_data); 
    D = reshape(data_compressed, size(data_compressed,1)*size(data_compressed,2), size(raw_data,3));
    [~,S,V] = svd(D,'econ');  
    ncoils_compressed = max(find(diag(S)/S(1) > 0.01)); % Keep >1% energy
    fprintf('Using %d compressed coils (out of %d original coils)\n', ...
            ncoils_compressed, size(raw_data,3));
    data_compressed = reshape(D * V(:,1:ncoils_compressed), ...
                   size(data_compressed,1), size(data_compressed,2), ncoils_compressed);

    % -- Prepare k-space data for reconstruction --

    % Upsample data
    upsample_factor = config.recon.upsample_factor;
    ns_new = (nsamples_per_petal - 1) * upsample_factor + 1;
    data_compressed_interp = zeros(ns_new, npetals, 3);
    x_old = 1:nsamples_per_petal;
    x_new = linspace(1, nsamples_per_petal, ns_new);
    
    for ndim = 1:3
        for p = 1:npetals
            y = data_compressed(:, p, ndim);     % original data
            data_compressed_interp(:, p, ndim) = interp1(x_old, y, x_new, config.recon.interpolation_method);
        end
    end

    % Separate data for te1 and te2 and skip and shift usampled k-space
    skip_start_te1 = config.recon.skip_start_te1;
    skip_end_te1 = config.recon.skip_end_te1;
    skip_start_te2 = config.recon.skip_start_te2;
    skip_end_te2 = config.recon.skip_end_te2;
    shift_te1 = config.recon.shift_te1;
    shift_te2 = config.recon.shift_te2;

    start_te1_upsampled = (start_te1 - 1) * upsample_factor + 1;
    start_te2_upsampled = (start_te2 - 1) * upsample_factor + 1;
    nsamples_per_echo_upsampled = (nsamples_per_petal/2 - 1)*upsample_factor + 1;

    data_te1_upsampled = data_compressed_interp(start_te1_upsampled + skip_start_te1 + shift_te1: ...
        start_te1_upsampled + nsamples_per_echo_upsampled - 1 - skip_end_te1 + shift_te1,:,:);
    data_te2_upsampled = data_compressed_interp(start_te2_upsampled + skip_start_te2 + shift_te2: ...
        start_te2_upsampled + nsamples_per_echo_upsampled - 1 - skip_end_te2 + shift_te2,:,:);
    
    % Downsample data after skipping and shifting
    data_te1 = data_te1_upsampled(1:upsample_factor:end, :, :);
    data_te2 = data_te2_upsampled(1:upsample_factor:end, :, :);
    
    % -- Reconstruction --

    % Check if valid density compensation file provided in config file
    if isfile(config.io.dcf_path)
        % Load density compensation weights
        disp('Density compensation weights file detected.');
        disp('Loading input density compensation weights...');
        load(config.io.dcf_path, 'dcf_weights_te1', 'dcf_weights_te2');

        % Estimate coil sensitivity map from first echo and recon first echo
        disp('Performing coilwise reconstruction for first echo...');
        coilwise_recon_te1 = coilwise_nufft_recon(data_te1, trajectory_te1, ...
            nx, ny, nz, dcf_weights_te1);
        disp('Estimating coil sensitivity map from first echo...');
        csm = estimate_csm(coilwise_recon_te1);
        disp('Perfoming coil combination for first echo...');
        img_recon_te1 = combine_coils(coilwise_recon_te1, csm);

        % Recon second echo using coil sensitivity map from first echo
        disp('Performing coilwise reconstruction for second echo...');
        coilwise_recon_te2 = coilwise_nufft_recon(data_te2, trajectory_te2, ...
            nx, ny, nz, dcf_weights_te2);
        disp('Perfoming coil combination for second echo...');
        img_recon_te2 = combine_coils(coilwise_recon_te2, csm);
    else
        % Generate density compensation weights during coilwise recon
        disp('No valid density compensation weights file detected');

        % Estimate coil sensitivity map from first echo and recon first echo
        disp('Generating density compensation weights and performing coilwise reconstruction for first echo...');
        coilwise_recon_te1 = coilwise_nufft_recon(data_te1, trajectory_te1, nx, ny,nz);
        disp('Estimating coil sensitivity map from first echo...');
        csm = estimate_csm(coilwise_recon_te1);
        img_recon_te1 = combine_coils(coilwise_recon_te1, csm);
    
        % Recon second echo using coil sensitivity map from first echo
        disp('Generating density compensation weights and performing coilwise reconstruction for second echo...');
        coilwise_recon_te2 = coilwise_nufft_recon(data_te2, trajectory_te2, nx, ny,nz);
        disp('Perfoming coil combination for second echo...');
        img_recon_te2 = combine_coils(coilwise_recon_te2, csm);
    end

    % -- Save outputs --
    % Create subdirectory for results
    out_dir = fullfile(config.io.out_path, 'UTERosette_bSSFP_dual_echo');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Save config, recon images, and trajectories as .mat files
    save(fullfile(out_dir, 'config.mat'), 'config');
    save(fullfile(out_dir, 'img_recon.mat'), 'img_recon_te1', 'img_recon_te2');
    save(fullfile(out_dir, 'trajectory.mat'), 'trajectory_te1', 'trajectory_te2');

    % Save initial image niftis to generate headers
    niftiwrite(img_recon_te1, fullfile(out_dir, 'img_recon_te1.nii')); 
    niftiwrite(img_recon_te2, fullfile(out_dir, 'img_recon_te2.nii'));
    
    % Read in nifti headers
    img_recon_te1_info = niftiinfo(fullfile(out_dir, 'img_recon_te1.nii'));
    img_recon_te2_info = niftiinfo(fullfile(out_dir, 'img_recon_te1.nii'));

    % Edit nifti header orientation and units information
    img_transform = eye(4); % define affine transformation for orientation
    img_transform(2,2) = -1; % reconstructed image is flipped along sagital plane
    img_recon_te1_info.Transform.T = img_transform;
    img_recon_te1_info.TransformName = 'Sform';
    img_recon_te1_info.SpaceUnits = 'Millimeter';
    img_recon_te2_info.Transform.T = img_transform;
    img_recon_te2_info.TransformName = 'Sform';
    img_recon_te2_info.SpaceUnits = 'Millimeter';

    % Resave (overwrite) niftis with edited headers
    niftiwrite(img_recon_te1, fullfile(out_dir, 'img_recon_te1.nii'), img_recon_te1_info);
    niftiwrite(img_recon_te2, fullfile(out_dir, 'img_recon_te2.nii'), img_recon_te2_info);

    total_processing_time = toc;
    fprintf('Total processing time: %.2f seconds\n', total_processing_time);

end