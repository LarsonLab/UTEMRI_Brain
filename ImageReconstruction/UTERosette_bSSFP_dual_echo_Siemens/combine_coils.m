function img_recon = combine_coils(x, csm_all)
%COMBINE_COILS Combines coil-wise images using provided sensitivity maps
%
% INPUTS:
%   x       : 4D array [nx, ny, nz, ncoils] of coil-wise reconstructed images
%   csm_all : 4D array [nx, ny, nz, ncoils] of coil sensitivity maps
%
% OUTPUT:
%   img_recon : 3D array [nx, ny, nz] of combined complex image

    [nx, ny, nz, ncoils] = size(x);
    img_recon = zeros(nx, ny, nz);

    for jj = 1:nz
        % Extract slice
        x_slice = squeeze(x(:,:,jj,:));       % [nx, ny, ncoils]
        csm_slice = squeeze(csm_all(:,:,jj,:)); % [nx, ny, ncoils]

        % Sum-of-squares normalization
        ll_sq = sum(csm_slice .* conj(csm_slice), 3);

        % Weighted combination
        img_recon(:,:,jj) = sum(conj(csm_slice) .* x_slice, 3) ./ ll_sq;
    end

end
