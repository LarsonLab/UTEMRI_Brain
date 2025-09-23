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

    % Startup/wind-down weighting
    weighted = linspace(0,0.95,44);
    weighted_cum = cumsum(weighted);
    
    % Initialize k-space arrays
    kx = zeros(190,378,432);
    ky = zeros(190,378,432);
    kz = zeros(190,378,432);
    
    % Build trajectory
    for ll = 1:190
        alpha = pi/189 * (ll - 1);
    
        for kk = 1:378
            phi = 2*pi/378 * (kk - 1);
    
            % Startup region
            for jj = 1:44
                w = weighted_cum(jj);
                ang = pi/(195*2) * w;
                kx(ll,kk,jj) = kmax*sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = kmax*sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = kmax*sin(ang)*sin(-pi/2+alpha);
            end
    
            % Main region
            for jj = 45:392
                ang = pi/(195*2) * (jj-23);
                kx(ll,kk,jj) = kmax*sin(ang)*cos(ang+phi)*cos(-pi/2+alpha);
                ky(ll,kk,jj) = kmax*sin(ang)*sin(ang+phi)*cos(-pi/2+alpha);
                kz(ll,kk,jj) = kmax*sin(ang)*sin(-pi/2+alpha);
            end
    
            % Wind-down region
            for jj = 393:432
                w = (195*2) - weighted_cum(437-jj);
                ang = pi/(195*2) * w;
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
    kx = pi / max(abs(kx(:))) * kx;
    ky = pi / max(abs(ky(:))) * ky;
    kz = pi / max(abs(kz(:))) * kz;

    % Extract first and second echo trajectories
    kx_te1 = kx(2:217,:,:);
    ky_te1 = ky(2:217,:,:);
    kz_te1 = kz(2:217,:,:);

    kx_te2 = kx(216:431,:,:);
    ky_te2 = ky(216:431,:,:);
    kz_te2 = kz(216:431,:,:);
    
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
