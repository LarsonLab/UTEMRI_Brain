function [trajectory_te1, trajectory_te2] = generate_dual_echo_rosette_trajectory(traj)
%GENERATE_DUAL_ECHO_ROSETTE_TRAJECTORY
% Generates a two 3D rosette k-space trajectories for each echo of dual echo UTE sequence scaled to ±pi.
%
% INPUTS:
%   traj   - Struct containing configurable trajectory parameters:
%       .npoints_petal      : Number of points per rosette petal
%       .npoints_ramp       : Number of points in non linear ramp up/down region of each petal
%       .ramp_ang           : Angle spanned by non linear ramp up/down region of each petal
%                             (NOTE: must be <=pi/2)
%       .ramp_weights       : [1 x npoints_ramp] array of weights to define ramp up/down scaling.
%                             (NOTE: should be monotonically increasing
%                             normalized to max(ramp_weights)=1
%       .npoints_alpha      : Number of points to sample alpha angle 
%                             (rotation of petals relative to the z axis through different polar angles)
%       .npoints_phi        : Number of points to sample phi angle 
%                             (rotation of petals in xy plane about the z axis)
%       .npoints_zero_start : Number of zero points in beginning of each petal
%       .start_te1          : Start point of first echo along each petal (setting a fraction will incrementally shift the trajectory)
%       .start_te2          : Start point of second echo along each petal (setting a fraction will incrementally shift the trajectory)
%
% OUTPUTS:
%   trajectory_te1 : trajectory for first echo. mx3 array of column vectors of normalized k-space
%                    coordinates in x, y, and z dimensions, respectively.
%   trajectory_te2 : trajectory for second echo. mx3 array of column vectors of normalized k-space
%                    coordinates in x, y, and z dimensions, respectively.

    % Calculate spacing of points in linear section of each petal
    linear_spacing = (pi - traj.ramp_ang*2)/(traj.npoints_petal ...
        - traj.npoints_ramp*2-traj.npoints_zero_start);
    
    % Scale ramp up/down weights
    ramp_weights = traj.ramp_weights * traj.ramp_ang/linear_spacing;
    
    % Set starting point and length of first and second echo along a petal
    start_te1 = traj.start_te1; % start of first echo
    start_te2 = traj.start_te2; % start of second echo
    npoints_te1 = traj.npoints_petal/2;
    npoints_te2 = traj.npoints_petal/2;
    
    % Calculate number of points in linear section of each petal
    npoints_linear = traj.npoints_petal - ...
        (traj.npoints_ramp*2 + traj.npoints_zero_start);
    
    % Calculate index for start of each region of petal
    idx_ramp_up_start = traj.npoints_zero_start + 1;
    idx_linear_start = idx_ramp_up_start + traj.npoints_ramp;
    idx_ramp_down_start = idx_linear_start + npoints_linear;
    
    % Initialize k-space arrays
    kx = zeros(traj.npoints_alpha,traj.npoints_phi,traj.npoints_petal);
    ky = zeros(traj.npoints_alpha,traj.npoints_phi,traj.npoints_petal);
    kz = zeros(traj.npoints_alpha,traj.npoints_phi,traj.npoints_petal);
    
    % Build trajectory
    for ll = 1:traj.npoints_alpha
        alpha = pi/(traj.npoints_alpha-1) * (ll - 1);
    
        for kk = 1:traj.npoints_phi
            phi = 2*pi/traj.npoints_phi * (kk - 1);
    
            % Initialize index counters for ramp weight array
            pp = 1;
            qq = 1;
    
            % Ramp-up region
            for jj = idx_ramp_up_start:idx_ramp_up_start + traj.npoints_ramp - 1
                weight = ramp_weights(pp);
                pp = pp + 1;
                
                ang = linear_spacing * weight;
                kx(ll,kk,jj) = sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = sin(ang)*sin(-pi/2+alpha);
            end
    
            % Linear region
            for jj = idx_linear_start:idx_linear_start + npoints_linear - 1
                weight = weight + 1;
                ang = linear_spacing * weight;
                
                kx(ll,kk,jj) = sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = sin(ang)*sin(-pi/2+alpha);
            end
    
            % Ramp-down region
            for jj = idx_ramp_down_start:idx_ramp_down_start + traj.npoints_ramp - 1
                if qq < length(ramp_weights)
                    weight = pi/linear_spacing - ramp_weights(end - qq);
                else
                    % repeat last point on trajectory
                    weight = pi/linear_spacing - ramp_weights(1);
                end
                qq = qq + 1;
    
                ang = linear_spacing * weight;
                kx(ll,kk,jj) = sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = sin(ang)*sin(-pi/2+alpha);
            end
        end
    end
    
    % Reorder dimensions
    kx = permute(kx,[3 2 1]);
    ky = permute(ky,[3 2 1]);
    kz = permute(kz,[3 2 1]);
    
    % Normalize to max = pi
    kx = pi / max(abs(kx(:))) * kx;
    ky = pi / max(abs(ky(:))) * ky;
    kz = pi / max(abs(kz(:))) * kz;
    
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
