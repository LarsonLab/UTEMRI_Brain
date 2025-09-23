clear;

data_file = '/working/larson3/Xin/VFA_MS_1/meas_MID03285_FID24049_fid_rosette_2mm_bSSFP_41k.dat';
write_flag = 0;  write_filename = 'testrecon';

%% Setup Dependencies
if (~exist('Gmri','file'))
    %set up IRT 
    currentDirectory = pwd;
    cd('/home/plarson/matlab/reconstruction/mirt/'); setup
    cd(currentDirectory)
end


%% Read data
[refmrprot, refmdh, inref] =rdMeas_dene(data_file);

%% Scan and recon parameters
sampling_interval=10; %in us
fov=0.24; %field of view in m
matrix_size=240; % number of samples along one spatial dimension

% reconstructed matrix size
nx=256;
ny=256;
nz=256;

ncoils = size(inref,3);
npetals = size(inref,2);
nsamples = size(inref,1);


%% Data modulations
% phase correction for chopping(?)
inref(:,2:2:end,:) = -inref(:,2:2:end,:);

% demodulate k-space data to account for off-resonance
frequency_offset = 300; %Hz
t=(0:nsamples-1) * sampling_interval/2;  % us
f_modulation = exp(1j*2*pi*frequency_offset*t'/1e6);
inref = inref .* repmat(f_modulation, [1, npetals, ncoils]);

%% Generate k-space trajectories

[trajectory_te1, trajectory_te2] = generate_dual_echo_rosette_trajectory(matrix_size, fov);


%% Reconstruct first echo

trajectory = trajectory_te1;
mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(trajectory,1); % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=(ir_mri_density_comp([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));

data_all=inref(3:218,:,:);
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[~,S,V] = svd(D,'econ');  
ncoils_compressed = max(find(diag(S)/S(1)>0.01)); %0.01
%%
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
ncoils_compressed
clear x;
for abc=1:ncoils_compressed
      lll=data(:,:,abc);
xcp=(Gm1'*(((lll(:)).*(w2)')));
	xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);

x(:,:,:,abc)=embed(xpcg,mask2(:,:,:));
abc
end

for jj=1:256
    px = angle(squeeze((x(:,:,jj,:))));

nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
ll=estimate_csm_walsh(nws_water_nuf);
ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
    
nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
img_recon(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
end

%%  Reconstruct second echo

trajectory = trajectory_te2;
mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(trajectory,1); % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=(ir_mri_density_comp([trajectory(:,1) trajectory(:,2) trajectory(:,3)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));

data=squeeze(inref(217:432,:,:));

D = reshape(data,size(data,1)*size(data,2),ncoils);
[~,S,V] = svd(D,'econ');  

data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
%%
ncoils_compressed
clear x;
for abc=1:ncoils_compressed
      lll=data(:,:,abc);
xcp=(Gm1'*(((lll(:)).*(w2)')));
	xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);

x(:,:,:,abc)=embed(xpcg,mask2(:,:,:));
abc
end

for jj=1:256 
    px = angle(squeeze((x(:,:,jj,:))));

nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
ll=estimate_csm_walsh(nws_water_nuf);
ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
    
nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
img_recon_TE2(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
end

%%
if write_flag
    save(write_filename,'img_recon_TE1', 'img_recon_TE2');
end