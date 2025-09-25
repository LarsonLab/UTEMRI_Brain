function [trajectory_te1, trajectory_te2] = generate_dual_echo_rosette_trajectory(matrix_size, fov)
%GENERATE_DUAL_ECHO_ROSETTE_TRAJECTORY
% Generates a two 3D rosette k-space trajectories for each echo of dual echo UTE sequence scaled to ±pi.
%
% INPUTS:
%   matrix_size : number of samples along one spatial dimension
%   fov         : field of view (in meters)
%
% OUTPUTS:
%   trajectory_te1 : trajectory for first echo. mx3 array of column vectors of normalized k-space
%                    coordinates in x, y, and z dimensions, respectively.
%   trajectory_te2 : trajectory for second echo. mx3 array of column vectors of normalized k-space
%                    coordinates in x, y, and z dimensions, respectively.

    % Calculate kmax
    kmax = matrix_size / fov / 2;

    % Set trajectory parameters
    npoints_petal = 432; % number of points per petal
    npoints_zero_start = 2; % number of zero points in beginning of each petal
    linear_spacing = pi/390; % spacing of points in linear section of each petal
    npoints_ramp = 40; % number of points in non linear ramp up/down of each petal
    npoints_alpha = 190; % number of points to sample alpha angle
    npoints_phi = 378; % number of points to sample phi angle
    ramp_weights = cumsum(linspace(0,1,npoints_ramp)); % defines the scaling of how trajectory is ramped up/down
    
    % Set starting point and length of first and second echo along a petal
    % (setting a fraction will incrementally shift the trajectory)
    start_te1 = 1; % start of first echo
    start_te2 = 214.65; % start of second echo
    npoints_te1 = npoints_petal/2;
    npoints_te2 = npoints_petal/2;
    
    % Calculate number of points in linear section of each petal
    npoints_linear = npoints_petal - (npoints_ramp*2 + npoints_zero_start); % number of points in linear section of each petal
    
    % Calculate index for start of each region of petal
    idx_ramp_up_start = npoints_zero_start + 1;
    idx_linear_start = idx_ramp_up_start + npoints_ramp;
    idx_ramp_down_start = idx_linear_start + npoints_linear;
    
    % Initialize k-space arrays
    kx = zeros(npoints_alpha,npoints_phi,npoints_petal);
    ky = zeros(npoints_alpha,npoints_phi,npoints_petal);
    kz = zeros(npoints_alpha,npoints_phi,npoints_petal);
    
    % Build trajectory
    for ll = 1:npoints_alpha
        alpha = pi/(npoints_alpha-1) * (ll - 1);
    
        for kk = 1:npoints_phi
            phi = 2*pi/npoints_phi * (kk - 1);
    
            % Initialize index counters for ramp weight array
            pp = 1;
            qq = 1;
    
            % Ramp-up region
            for jj = idx_ramp_up_start:idx_ramp_up_start + npoints_ramp - 1
                weight = ramp_weights(pp);
                pp = pp + 1;
                
                ang = linear_spacing * weight;
                kx(ll,kk,jj) = kmax*sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = kmax*sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = kmax*sin(ang)*sin(-pi/2+alpha);
            end
    
            % Linear region
            for jj = idx_linear_start:idx_linear_start + npoints_linear - 1
                weight = weight + 1;
                ang = linear_spacing * weight;
                
                kx(ll,kk,jj) = kmax*sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = kmax*sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = kmax*sin(ang)*sin(-pi/2+alpha);
            end
    
            % Ramp-down region
            for jj = idx_ramp_down_start:idx_ramp_down_start + npoints_ramp - 1
                if qq < length(ramp_weights)
                    weight = pi/linear_spacing - ramp_weights(end - qq);
                else
                    % repeat last point on trajectory
                    weight = pi/linear_spacing - ramp_weights(1);
                end
                qq = qq + 1;
    
                ang = pi/(195*2) * weight;
                kx(ll,kk,jj) = kmax*sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = kmax*sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = kmax*sin(ang)*sin(-pi/2+alpha);
            end
        end
    end
    
    % Reorder dimensions
    kx = permute(kx,[3 2 1]);
    ky = permute(ky,[3 2 1]);
    kz = permute(kz,[3 2 1]);
    
    % Normalize to max = pi
    % kx = pi / max(abs(kx(:))) * kx;
    % ky = pi / max(abs(ky(:))) * ky;
    % kz = pi / max(abs(kz(:))) * kz;
    
    % Extract first and second echo trajectories (including shifting)
    shift_frac_te1 = start_te1 - floor(start_te1);
    kx_te1 = (1-shift_frac_te1)*kx(floor(start_te1):floor(start_te1)+npoints_te1-1,:,:) + ...
        shift_frac_te1*kx(floor(start_te1)+1:floor(start_te1)+npoints_te1,:,:);
    ky_te1 = (1-shift_frac_te1)*ky(floor(start_te1):floor(start_te1)+npoints_te1-1,:,:) + ...
        shift_frac_te1*ky(floor(start_te1)+1:floor(start_te1)+npoints_te1,:,:);
    kz_te1 = (1-shift_frac_te1)*kz(floor(start_te1):floor(start_te1)+npoints_te1-1,:,:) + ...
        shift_frac_te1*kz(floor(start_te1)+1:floor(start_te1)+npoints_te1,:,:);
    
    shift_frac_te2 = start_te2 - floor(start_te2);
    kx_te2 = (1-shift_frac_te2)*kx(floor(start_te2):floor(start_te2)+npoints_te2-1,:,:) + ...
        shift_frac_te2*kx(floor(start_te2)+1:floor(start_te2)+npoints_te2,:,:);
    ky_te2 = (1-shift_frac_te2)*ky(floor(start_te2):floor(start_te2)+npoints_te2-1,:,:) + ...
        shift_frac_te2*ky(floor(start_te2)+1:floor(start_te2)+npoints_te2,:,:);
    kz_te2 = (1-shift_frac_te2)*kz(floor(start_te2):floor(start_te2)+npoints_te2-1,:,:) + ...
        shift_frac_te2*kz(floor(start_te2)+1:floor(start_te2)+npoints_te2,:,:);

    % Normalize to max = pi
    kx_te1 = pi / max(abs(kx_te1(:))) * kx_te1;
    ky_te1 = pi / max(abs(ky_te1(:))) * ky_te1;
    kz_te1 = pi / max(abs(kz_te1(:))) * kz_te1;
    
    kx_te2 = pi / max(abs(kx_te2(:))) * kx_te2;
    ky_te2 = pi / max(abs(ky_te2(:))) * ky_te2;
    kz_te2 = pi / max(abs(kz_te2(:))) * kz_te2;
        
    % Store trajectories in 4D arrays
    trajectory_te1 = zeros(length(kx_te1(:)), 3);
    trajectory_te2 = zeros(length(kx_te2(:)), 3);

    trajectory_te1(:,1) = kx_te1(:);
    trajectory_te1(:,2) = ky_te1(:);
    trajectory_te1(:,3) = kz_te1(:);

    trajectory_te2(:,1) = kx_te2(:);
    trajectory_te2(:,2) = ky_te2(:);
    trajectory_te2(:,3) = kz_te2(:);
end
