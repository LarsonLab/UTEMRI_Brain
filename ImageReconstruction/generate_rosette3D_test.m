function traj = generate_rosette3D_test(varargin)
% GENERATE_ROSETTE3D Generate 3D rosette k-space trajectories (Shen et al. MRM 2024)
%
% USAGE
%  traj = generate_rosette3D()                      % use paper defaults
%  traj = generate_rosette3D('Np',71820,'kmax',500) % override params
%
% OUTPUT (struct)
%  traj.kx, traj.ky, traj.kz    - arrays of size [Ns x Np] (k in m^-1)
%  traj.kx_vec, traj.ky_vec     - concatenated column vectors [Ns*Np x 1]
%  traj.t                        - time vector for a single petal (s)
%  traj.Gx, traj.Gy, traj.Gz    - gradient waveforms [Ns x Np] (T/m) (optional)
%  traj.params                   - used parameters
%
% DEFAULTS (matching Shen et al. 2024)
%  omega1 = omega2 = 1611;      % rad/s
%  kmax   = 500;                % m^-1
%  Ns     = 432;                % samples per petal
%  readout_dur = 2.16e-3;       % seconds (paper readout duration)
%  Np     = 2000;               % number of petals (use 71820 for full-coverage)
%  phi sampled uniform in [-pi/2,pi/2]
%  beta sampled uniform in [0,2*pi)
%
% Note: If you want to exactly reproduce the paper's full dataset set Np=71820.
%
% Relationship used for gradients:
%   k(t) = (gamma/2/pi) * integral(G dt)  => G = (2*pi/gamma) * dk/dt
%   Using gamma/(2*pi) = 42.576e6 Hz/T => G = dk/dt / (42.576e6)
%
% Example:
%   traj = generate_rosette3D('Np',1000,'plotOn',true);
%
% Author: ChatGPT (adapted to Shen et al. 2024)
% Date: 2025-09-18

% -------------------------
% Parse inputs / defaults
% -------------------------
p = inputParser;
addParameter(p,'omega1',1611,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'omega2',1611,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'kmax',500,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'Ns',432,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'readout_dur',2.16e-3,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'Np',2000,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'phi_vals',[],@(x)isnumeric(x));
addParameter(p,'beta_vals',[],@(x)isnumeric(x));
addParameter(p,'seed',0,@(x)isnumeric(x));
addParameter(p,'computeGradients',true,@islogical);
addParameter(p,'plotOn',false,@islogical);
parse(p,varargin{:});
pr = p.Results;

rng(pr.seed);

% time vector for single petal
Ns = pr.Ns;
t = linspace(0, pr.readout_dur, Ns)';    % s, column vector
omega1 = pr.omega1;
omega2 = pr.omega2;
kmax = pr.kmax;
Np = pr.Np;

% phi and beta sampling
if isempty(pr.phi_vals)
    % uniform sampling in [-pi/2, +pi/2]
    pr.phi_vals = linspace(-pi/2, +pi/2, Np);
end
if isempty(pr.beta_vals)
    % uniform sample of beta in [0, 2pi)
    pr.beta_vals = 2*pi * (0:(Np-1))'/Np;
end
if numel(pr.phi_vals) ~= Np
    error('phi_vals must have length Np (or be empty).');
end
if numel(pr.beta_vals) ~= Np
    error('beta_vals must have length Np (or be empty).');
end

% Preallocate per-petal arrays
kx = zeros(Ns, Np);
ky = zeros(Ns, Np);
kz = zeros(Ns, Np);

% Compute sin(omega1 t) once
s1 = sin(omega1 * t);        % Ns x 1

% Generate trajectories (vectorized across petals)
% Equation (paper):
% kxy(t) = (kmax*cos(phi)) * sin(omega1 t) * exp(i*(omega2*t + beta))
% kz(t)  = kmax*sin(phi) * sin(omega1 t)
% We'll evaluate real and imag parts of kxy -> kx, ky
for pidx = 1:Np
    phi = pr.phi_vals(pidx);
    beta = pr.beta_vals(pidx);
    amp_xy = kmax * cos(phi);
    ang = omega2 * t + beta;
    kxy = amp_xy .* s1 .* exp(1i * ang);   % Ns x 1 complex
    kx(:,pidx) = real(kxy);
    ky(:,pidx) = imag(kxy);
    kz(:,pidx) = kmax * sin(phi) .* s1;
end

% Concatenate into single vectors (trajectory streaming order: petals in index order)
kx_vec = kx(:);
ky_vec = ky(:);
kz_vec = kz(:);

% Optionally compute gradient waveforms (finite-difference)
Gx = []; Gy = []; Gz = [];
if pr.computeGradients
    % gamma/(2pi) = 42.576e6 Hz/T  => G = dk/dt / (42.576e6)
    gamma_over_2pi = 42.576e6; % Hz/T
    % compute dk/dt along time dimension for each petal
    dt = mean(diff(t));
    dkx_dt = gradient(kx, dt);   % Ns x Np
    dky_dt = gradient(ky, dt);
    dkz_dt = gradient(kz, dt);
    Gx = dkx_dt / gamma_over_2pi;    % T/m
    Gy = dky_dt / gamma_over_2pi;
    Gz = dkz_dt / gamma_over_2pi;
end

% Pack outputs
traj.kx = kx;
traj.ky = ky;
traj.kz = kz;
traj.kx_vec = kx_vec;
traj.ky_vec = ky_vec;
traj.kz_vec = kz_vec;
traj.t = t;
traj.Gx = Gx;
traj.Gy = Gy;
traj.Gz = Gz;
traj.params = pr;

% Optional plotting
if pr.plotOn
    figure('Name','3D Rosette k-space (subset)');
    % plot subset of petals for clarity
    plotN = min(200, Np);
    hold on;
    for ii = 1:plotN
        plot3(kx(:,ii), ky(:,ii), kz(:,ii), '-');
    end
    axis equal; grid on;
    xlabel('k_x (m^{-1})'); ylabel('k_y (m^{-1})'); zlabel('k_z (m^{-1})');
    title(sprintf('3D Rosette k-space (first %d petals)', plotN));
    view(3);
    hold off;
    
    % show center crossing for a single petal
    figure('Name','Single petal kx/ky/kz vs time');
    subplot(3,1,1); plot(t*1e3, kx(:,1)); xlabel('t (ms)'); ylabel('k_x (m^{-1})');
    subplot(3,1,2); plot(t*1e3, ky(:,1)); xlabel('t (ms)'); ylabel('k_y (m^{-1})');
    subplot(3,1,3); plot(t*1e3, kz(:,1)); xlabel('t (ms)'); ylabel('k_z (m^{-1})');
end

end
