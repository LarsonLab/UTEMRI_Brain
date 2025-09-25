clear;

%% Select input and output configurations
write_flag = 0; % 1 - save outputs, 0 - don't save output
load_trajectory_flag = 0; % 1 - load existing trajectory, 0 - generate trajectory

twix_path = ''; % path to twix file
trajectory_path = ''; % path to .mat trajectory file containing trajectory_te1 and trajectory_te2

out_path = ''; % path to save output files to

%% Set up Dependencies
mirt_path = '/home/plarson/matlab/reconstruction/mirt/'; % path to mirt reconstruction toolbox

if (~exist('Gmri','file'))
    %set up IRT 
    currentDirectory = pwd;
    cd(mirt_path);
    setup
    cd(currentDirectory)
end

%% Read in raw data
[refmrprot, refmdh, inref] =rdMeas_dene(twix_path);

%% Set scan and recon parameters
recon_parameters.sampling_interval=10; %in us
recon_parameters.fov=0.24; %field of view in m
recon_parameters.matrix_size=240; % number of samples along one spatial dimension
recon_parameters.npoints_skip_te1 = 2; % number of data points to skip in first echo
recon_parameters.npoints_skip_te2 = 0; % number of data points to skip in second echo

% reconstructed matrix size
nx=256;
ny=256;
nz=256;

recon_parameters.ncoils = size(inref,3);
recon_parameters.npetals = size(inref,2);
recon_parameters.nsamples_per_petal = size(inref,1);

% calculate starting data point of each echo (along each petal in rosette)
recon_parameters.start_te1 = recon_parameters.npoints_skip_te1 + 1;
recon_parameters.start_te2 = recon_parameters.nsamples_per_petal/2 + recon_parameters.npoints_skip_te2 + 1;

%% Data modulations
% phase correction for chopping(?)
inref(:,2:2:end,:) = -inref(:,2:2:end,:);

% demodulate k-space data to account for off-resonance
recon_parameters.frequency_offset = 300; %Hz
t=(0:recon_parameters.nsamples_per_petal-1) * recon_parameters.sampling_interval/2;  % us
recon_parameters.f_modulation = exp(1j*2*pi*recon_parameters.frequency_offset*t'/1e6);
inref = inref .* repmat(recon_parameters.f_modulation, [1, recon_parameters.npetals, recon_parameters.ncoils]);

%% Get k-space trajectories

if load_trajectory_flag
    load(trajectory_path, 'trajectory_te1', 'trajectory_te2');
    recon_parameters.trajectory_te1 = trajectory_te1;
    recon_parameters.trajectory_te2 = trajectory_te2;
else
    [recon_parameters.trajectory_te1, recon_parameters.trajectory_te2] = ...
        generate_dual_echo_rosette_trajectory(recon_parameters.matrix_size, recon_parameters.fov);
end

%% Reconstruct first echo
trajectory = recon_parameters.trajectory_te1;
mask = true([nx ny nz]);
nufft = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi), mask, 'fov', 256, 'basis', f.basis, 'nufft', nufft);
beta = 2^-21 * size(trajectory,1); % good for quadratic 'rect'
R = Reg1(mask, 'beta', beta);
kdens=(ir_mri_density_comp([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));

data=squeeze(inref(recon_parameters.start_te1:recon_parameters.start_te1 + recon_parameters.nsamples_per_petal/2 - 1,:,:));
D = reshape(data,size(data,1)*size(data,2),recon_parameters.ncoils);
[~,S,V] = svd(D,'econ');  
ncoils_compressed = max(find(diag(S)/S(1)>0.01)); %0.01
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);

x = zeros(size(mask), ncoils_compressed);
for ncoil=1:ncoils_compressed
    lll=data(:,:,ncoil);
    xcp=(Gm1'*(((lll(:)).*(w2)')));
    xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);
    x(:,:,:,ncoil)=embed(xpcg,mask(:,:,:));
    ncoil
end

img_recon = zeros(nx, ny, nz);
for jj=1:nz
    px = angle(squeeze((x(:,:,jj,:))));

    nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
    ll=estimate_csm_walsh(nws_water_nuf);
    ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
        
    nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
    img_recon(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
end

img_recon_te1 = flip(img_recon,1);

%%  Reconstruct second echo
trajectory = recon_parameters.trajectory_te2;
mask = true([nx ny nz]);
nufft = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi), mask, 'fov', 256, 'basis', f.basis, 'nufft', nufft);
beta = 2^-21 * size(trajectory,1); % good for quadratic 'rect'
R = Reg1(mask, 'beta', beta);
kdens=(ir_mri_density_comp([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));

data=squeeze(inref(recon_parameters.start_te2:recon_parameters.start_te2 + recon_parameters.nsamples_per_petal/2 - 1,:,:));
D = reshape(data,size(data,1)*size(data,2),recon_parameters.ncoils);
[~,S,V] = svd(D,'econ');  
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);

x = zeros(size(mask), ncoils_compressed);
for ncoil=1:ncoils_compressed
    lll=data(:,:,ncoil);
    xcp=(Gm1'*(((lll(:)).*(w2)')));
    xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);
    x(:,:,:,ncoil)=embed(xpcg,mask(:,:,:));
    ncoil
end

img_recon = zeros(nx, ny, nz);
for jj=1:nz
    px = angle(squeeze((x(:,:,jj,:))));

    nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
    ll=estimate_csm_walsh(nws_water_nuf);
    ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
        
    nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
    img_recon(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
end

img_recon_te2 = flip(img_recon,1);

%% Save outputs
[~, twix_filename, file_ext] = fileparts(twix_path);
recon_parameters.twix_filename = [twix_filename, file_ext];

if write_flag
    save([out_path 'img_recon.mat'],'img_recon_te1', 'img_recon_te2');
    niftiwrite(img_recon_te1, [out_path 'img_recon_te1.nii']);
    niftiwrite(img_recon_te2, [out_path 'img_recon_te2.nii']);
    save([out_path 'recon_parameters.mat'], 'recon_parameters');
end