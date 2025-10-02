function config = base_config()
    %BASE_CONFIG Returns a struct with default fields and parameters for
    %UTE Rosette dual echo reconstruction

    % Input/output fields
    config.io.twix_path = ''; % path to raw Siemens twix file
    config.io.out_path = ''; % path to save outputs to

    % Reconstruction parameters
    config.recon.sampling_interval=10; %in us
    config.recon.fov=0.24; %field of view in m
    config.recon.matrix_size=240; % matrix size used for trajectory generation
    config.recon.npoints_skip_te1 = 2; % number of data points to skip in first echo
    config.recon.npoints_skip_te2 = 0; % number of data points to skip in second echo
    config.recon.output_size.nx = 256; % reconstructed matrix x size
    config.recon.output_size.ny = 256; % reconstructed matrix y size
    config.recon.output_size.nz = 256; % reconstructed matrix z size

    % Trajectory parameters
    config.recon.traj.npoints_petal = 432; % number of points per rosette petal
    config.recon.traj.npoints_ramp = 40; % number of points in non linear ramp up/down region of each petal
    config.recon.traj.ramp_ang = 20*pi/390; % angle spanned by non linear ramp up/down region of each petal
    ramp_weights = cumsum(linspace(0.05,0.95,config.recon.traj.npoints_ramp));
    config.recon.traj.ramp_weights = ramp_weights/ramp_weights(end); % [1 x npoints_ramp] 
                                                                     % array of weights to define ramp up/down scaling.
    config.recon.traj.npoints_alpha = 190; % number of points to sample alpha angle 
                                           % (rotation of petals relative to the z axis through different polar angles)
    config.recon.traj.npoints_phi = 378; % number of points to sample phi angle (rotation of petals in xy plane about the z axis)
    config.recon.traj.npoints_zero_start = 2; % number of zero points in beginning of each petal
    config.recon.traj.start_te1 = 1; % start point of first echo
    config.recon.traj.start_te2 = 214.65; % start point of second echo
end