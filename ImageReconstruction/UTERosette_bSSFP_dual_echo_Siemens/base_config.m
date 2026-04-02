function config = base_config()
    %BASE_CONFIG Returns a struct with default fields and parameters for
    %UTE Rosette dual echo reconstruction

    % NOTES:
    %   trajectory_te1  : trajectory for first echo. mx3 array of column vectors of normalized k-space
    %                     coordinates in x, y, and z dimensions, respectively.
    %   trajectory_te2  : trajectory for second echo. mx3 array of column vectors of normalized k-space
    %                     coordinates in x, y, and z dimensions, respectively.
    %   dcf_weights_te1 : density compensation function weights for first
    %                     echo
    %   dcf_weights_te2 : density compensation function weights for second
    %                     echo


    % Input/output fields
    config.io.twix_path = ''; % path to raw Siemens twix file
    config.io.out_path = ''; % path to save outputs to
    config.io.trajectory_path = ''; % path to .mat file containing trajectories (trajectory_te1, trajectory_te2)
    config.io.dcf_path = ''; % path to .mat file containing dcf weights (dcf_weights_te1, dcf_weights_te2)

    % Reconstruction parameters
    config.recon.sampling_interval=10; %in us
    config.recon.fov=0.24; %field of view in m
    config.recon.matrix_size=240; % matrix size used for trajectory generation
    config.recon.output_size.nx = 256; % reconstructed matrix x size
    config.recon.output_size.ny = 256; % reconstructed matrix y size
    config.recon.output_size.nz = 256; % reconstructed matrix z size
    config.recon.frequency_offset = 300; % frequency offset used to demodulate kspace data (Hz)
    config.recon.upsample_factor = 10; % factor by which to increase kspace data points through interpolation
    config.recon.interpolation_method = 'spline'; % method by which to interpolate kspace data for upsampling
    config.recon.skip_start_te1 = 0; % number of upsampled te1 kspace points to skip at the start of each petal
    config.recon.skip_end_te1 = 0; % number of upsampeld te1 kspace points to skip at the end of each petal
    config.recon.skip_start_te2 = 0; % number of upsampled te2 kspace points to skip at the start of each petal
    config.recon.skip_end_te2 = 0; % number of upsampeld te2 kspace points to skip at the end of each petal
    config.recon.shift_te1 = 2*config.recon.upsample_factor; % number of upsampled te1 kspace points to shift along each petal
    config.recon.shift_te2 = 0; % number of upsampled te2 kspace points to shift along each petal

    % Trajectory parameters
    config.recon.traj.npoints_petal = 432; % number of points per rosette petal
    config.recon.traj.npoints_ramp = 40; % number of points in non linear ramp up/down region of each petal
    config.recon.traj.ramp_ang = 20*pi/390; % angle spanned by non linear ramp up/down region of each petal
    ramp_weights = cumsum(linspace(0.05,0.95,config.recon.traj.npoints_ramp));
    config.recon.traj.ramp_weights = ramp_weights/ramp_weights(end); % [1 x npoints_ramp] array of weights to define ramp up/down scaling.
    config.recon.traj.npoints_alpha = 190; % number of points to sample alpha angle (rotation of petals relative to the z axis through different polar angles)
    config.recon.traj.npoints_phi = 378; % number of points to sample phi angle (rotation of petals in xy plane about the z axis)
    config.recon.traj.npoints_zero_start = 2; % number of zero points in beginning of each petal
    config.recon.traj.start_te1 = 1; % start point of first echo
    config.recon.traj.start_te2 = 214.65; % start point of second echo
    config.recon.traj.skip_start_te1 = 0; % number of points in each petal to skip at start of te1 trajectory
    config.recon.traj.skip_end_te1 = 0; % number of points in each petal to skip at end of te1 trajectory
    config.recon.traj.skip_start_te2 = 0; % number of points in each petal to skip at start of te2 trajectory
    config.recon.traj.skip_end_te2 = 0; % number of points in each petal to skip at end of te2 trajectory
end