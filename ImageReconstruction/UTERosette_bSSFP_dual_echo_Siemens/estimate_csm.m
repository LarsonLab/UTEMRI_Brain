function csm_all = estimate_csm(x)
%ESTIMATE_CSM Estimates coil sensitivity maps for all slices
%
% INPUT:
%   x       : 4D array [nx, ny, nz, ncoils] of complex coil images
%
% OUTPUT:
%   csm_all : 4D array [nx, ny, nz, ncoils] of complex coil sensitivity maps

    [nx, ny, nz, ncoils] = size(x);
    csm_all = zeros(nx, ny, nz, ncoils);

    tic;
    parfor jj = 1:nz
        % Extract slice
        x_slice = squeeze(x(:,:,jj,:));   % [nx, ny, ncoils]

        % Remove phase
        px = angle(x_slice);
        x_phase_corrected = x_slice .* exp(-1i * px);

        % Estimate coil sensitivities (Walsh method)
        ll = estimate_csm_walsh(x_phase_corrected);
        ll(ll < eps) = 1;

        % Store
        csm_all(:,:,jj,:) = ll;
    end
    csm_estimation_time = toc;
    fprintf('Total time for coil sensitivity map estimation: %.2f seconds\n', csm_estimation_time);

end
