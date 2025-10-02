function img_recon = reconstruct_coil_compressed(data, trajectory, nx, ny, nz)
%RECONSTRUCT_COIL_COMPRESSED  Reconstructs an MRI image from non-Cartesian k-space
%   using NUFFT-based reconstruction with coil compression.
%
%   INPUTS:
%       data             - Raw k-space data array, typically size 
%       recon_parameters - Struct with reconstruction parameters:
%       trajectory       - non-Cartesian k-space trajectory. [npoints x 3] array of column vectors of k-space
%                          coordinates in x, y, and z dimensions, respectively.
%       nx, ny, nz       - Output image matrix size
%
%   OUTPUT:
%       img_recon        - Final reconstructed, coil-combined image [nx x ny x nz]
%
%   Dependencies:
%       - MIRT toolbox
%       - Walsh coil sensitivity estimation (estimate_csm_walsh)

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
    
    % Coil compression via SVD
    data_compressed = squeeze(data); 
    D = reshape(data_compressed, size(data_compressed,1)*size(data_compressed,2), size(data,3));
    [~,S,V] = svd(D,'econ');  
    ncoils_compressed = max(find(diag(S)/S(1) > 0.01)); % Keep >1% energy
    fprintf('Using %d compressed coils (out of %d original coils)\n', ...
            ncoils_compressed, size(data,3));
    data_compressed = reshape(D * V(:,1:ncoils_compressed), ...
                   size(data_compressed,1), size(data_compressed,2), ncoils_compressed);
    
    % Coil-wise image reconstruction
    x = zeros(nx, ny, nz, ncoils_compressed);   
    for ncoil = 1:ncoils_compressed
        lll = data_compressed(:,:,ncoil);
        xcp = (Gm1' * ((lll(:)) .* (w)'));
        xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);
        x(:,:,:,ncoil) = embed(xpcg, mask(:,:,:));
        ncoil
    end
    
    % Coil combination using Walsh method
    img_recon = zeros(nx, ny, nz);
    for jj = 1:nz
        px = angle(squeeze(x(:,:,jj,:)));
        
        % Remove phase
        nws_water_nuf = squeeze(x(:,:,jj,:)) .* exp(-1i * px);
        
        % Estimate coil sensitivities
        ll = estimate_csm_walsh(nws_water_nuf);
        ll_sq = sum(ll .* conj(ll), 3); 
        ll(ll < eps) = 1;
        
        % Apply coil combination
        nws_water_nuf = squeeze(x(:,:,jj,:)) .* exp(-1i * px);
        img_recon(:,:,jj) = sum(conj(ll) .* nws_water_nuf ,3) ./ ll_sq;
    end
    
    % Flip reconstructed image
    img_recon = flip(img_recon,1);

end
