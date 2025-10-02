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

end