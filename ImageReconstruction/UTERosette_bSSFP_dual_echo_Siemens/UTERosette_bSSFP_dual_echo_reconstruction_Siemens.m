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
%       .recon.output_size: Reconstructed matrix size
%                           (fields: nx, ny, nz)
%       .recon.matrix_size: Matrix size used for trajectory generation
%       .recon.fov        : Field of view (mm)
%       .recon.npoints_skip_te1 : # of samples to skip at start of TE1
%       .recon.npoints_skip_te2 : # of samples to skip at start of TE2
%       .recon.sampling_interval: Sampling interval (µs)
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

    % -- Read in raw data --
    [~, ~, raw_data] =rdMeas_dene(config.io.twix_path);
    
    % -- Setup some recon parameters --
    % reconstructed matrix size
    nx=config.recon.output_size.nx;
    ny=config.recon.output_size.ny;
    nz=config.recon.output_size.nz;
    
    ncoils = size(raw_data,3);
    npetals = size(raw_data,2);
    nsamples_per_petal = size(raw_data,1);
    
    % calculate starting data point of each echo (along each petal in rosette)
    start_te1 = config.recon.npoints_skip_te1 + 1;
    start_te2 = nsamples_per_petal/2 + config.recon.npoints_skip_te2 + 1;
    
    % -- Data modulations --
    % phase correction for chopping(?)
    raw_data(:,2:2:end,:) = -raw_data(:,2:2:end,:);
    
    % demodulate k-space data to account for off-resonance
    recon_parameters.frequency_offset = 300; %Hz
    t=(0:nsamples_per_petal-1) * config.recon.sampling_interval/2;  % us
    recon_parameters.f_modulation = exp(1j*2*pi*recon_parameters.frequency_offset*t'/1e6);
    raw_data = raw_data .* repmat(recon_parameters.f_modulation, [1, npetals, ncoils]);
    
    % -- Create k-space trajectories --
    [trajectory_te1, trajectory_te2] = ...
        generate_dual_echo_rosette_trajectory(config.recon.traj);
    
    % -- Reconstruction --
    % first echo
    data=squeeze(raw_data(start_te1:start_te1 + nsamples_per_petal/2 - 1,:,:));
    img_recon_te1 = reconstruct_coil_compressed(data, trajectory_te1, nx, ny,nz);
    
    % second echo
    data=squeeze(raw_data(start_te2:start_te2 + nsamples_per_petal/2 - 1,:,:));
    img_recon_te2 = reconstruct_coil_compressed(data, trajectory_te2, nx, ny,nz);


    % -- Save outputs --
    % Create subdirectory for results
    out_dir = fullfile(config.io.out_path, 'UTERosette_bSSFP_dual_echo');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Save config, recon images, and trajectories
    save(fullfile(out_dir, 'config.mat'), 'config');
    save(fullfile(out_dir, 'img_recon.mat'), 'img_recon_te1', 'img_recon_te2');
    save(fullfile(out_dir, 'trajectory.mat'), 'trajectory_te1', 'trajectory_te2');
    niftiwrite(img_recon_te1, fullfile(out_dir, 'img_recon_te1.nii'));
    niftiwrite(img_recon_te2, fullfile(out_dir, 'img_recon_te2.nii'));
   
end