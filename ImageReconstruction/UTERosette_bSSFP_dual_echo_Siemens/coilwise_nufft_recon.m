function coilwise_recon = coilwise_nufft_recon(data, trajectory, nx, ny, nz)
%COILWISE_NUFFT_RECON Reconstructs an MRI image from non-Cartesian k-space
%   using NUFFT-based reconstruction.
%
%   INPUTS:
%       data             - Raw k-space data array, typically size 
%       trajectory       - non-Cartesian k-space trajectory. [npoints x 3] array of column vectors of k-space
%                          coordinates in x, y, and z dimensions, respectively.
%       nx, ny, nz       - Output image matrix size
%
%   OUTPUT:
%       coilwise_recon   - Reconstructed image from each coil [nx x ny x nz x ncoils]
%
%   Dependencies:
%       - MIRT toolbox

    % NUFFT and regularization setup
    mask = true([nx ny nz]);    
    nufft = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2], ...
             'table', 2^10, 'minmax:kb'};  
    f.basis = {'rect'};
    Gm1 = Gmri(trajectory/(2*pi), mask, 'fov', 256, ...
               'basis', f.basis, 'nufft', nufft);  
    beta = 2^-21 * size(trajectory,1); % Regularization weight
    R = Reg1(mask, 'beta', beta);
    
    % Density compensation
    kdens = ir_mri_density_comp(trajectory/(2*pi), 'pipe', ...
                                'G', Gm1.arg.Gnufft, 'arg_pipe', {'fov',256})';
    w = kdens / max(kdens(:));
    
    % Coil-wise image reconstruction
    ncoils = size(data, 3);
    coilwise_recon = zeros(nx, ny, nz, ncoils);   
    for ncoil = 1:ncoils
        lll = data(:,:,ncoil);
        xcp = (Gm1' * ((lll(:)) .* (w)'));
        xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);
        coilwise_recon(:,:,:,ncoil) = embed(xpcg, mask(:,:,:));
        ncoil
    end

end