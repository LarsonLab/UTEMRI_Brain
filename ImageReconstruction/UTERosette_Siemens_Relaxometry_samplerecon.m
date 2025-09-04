
%%
data_directory = '/working/larson3/reVIVE/revive_5-31-24/';
data_file = 'meas_MID03131_FID03159_fid_rosette_2mm2x_TE40.dat';
% '/working/larson2/xin/Revive_03_22_24/meas_MID00125_FID23366_fid_rosette_2mm2x_TE40.dat'

write_flag = 0;  write_filename = 'testrecon';

%% Setup Dependencies
if (~exist('Gmri','file'))
    %set up IRT 
    currentDirectory = pwd;
    cd('/home/plarson/matlab/reconstruction/mirt/'); setup
    cd(currentDirectory)
end

%% Read data

[refmrprot, refmdh, inref] =rdMeas_dene([data_directory data_file]);

%% scan and reconstruction parameters
gamma=42.576; %kHZ/mT
TR=2500; %us
sampling_interval=10; %us
FOV=0.24; %m
T=1.9; %ms

matrix_size=120; %120x120 isotropic for trajectory

ncoils  =   size(inref,3);

% nufft matrix size
nx=128;
ny=128;
nz=128;


%% k-space trajectory
Kmax=matrix_size/FOV/2*2*pi;
K_interval=Kmax*2/matrix_size;
n=1:1:60;
kr=n*K_interval;
Gr=kr*2*pi/(gamma*T);
dead_time=0.2;%ms
Gx_initial=2*kr/gamma/dead_time;
Gy_initial=-Gr;
time_pre=0.2;%ms
weighted=[0.01,0.03,0.05,0.08,0.12,0.2,0.37,0.55,0.72,0.87];
weighted_cum=cumsum(weighted);
for ll=1:95
    alpha=pi/189*2*(ll-1);
    for kk=1:190
        kx(ll,kk,1)=0;
        ky(ll,kk,1)=0;
        kz(ll,kk,1)=0;
        Gx(ll,kk,1)=0;
        Gy(ll,kk,1)=0;
        Gz(ll,kk,1)=0;
        kx(ll,kk,2)=0;
        ky(ll,kk,2)=0;
        kz(ll,kk,2)=0;
        Gx(ll,kk,2)=0;
        Gy(ll,kk,2)=0;
        Gz(ll,kk,2)=0;
        phi=2*pi/190*(kk-1);
        phi_next=2*pi/190*kk;
        for jj=1:10
            kx(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*cos(pi/195*weighted_cum(jj)+phi)*cos(-pi/2+alpha);
            ky(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*sin(pi/195*weighted_cum(jj)+phi)*cos(-pi/2+alpha);
            kz(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*sin(-pi/2+alpha);
        end
        for ii=4:191
            kx(ll,kk,ii-3+11+1)=Kmax*sin(pi/195*ii)*cos(pi/195*ii+phi)*cos(-pi/2+alpha);
            ky(ll,kk,ii-3+11+1)=Kmax*sin(pi/195*ii)*sin(pi/195*ii+phi)*cos(-pi/2+alpha);
            kz(ll,kk,ii-3+11+1)=Kmax*sin(pi/195*ii)*sin(-pi/2+alpha);
        end
        for jj=1:10
            kx(ll,kk,jj+199+1)=Kmax*sin(pi/195*(195-weighted_cum(11-jj)))*cos(pi/195*(195-weighted_cum(11-jj))+phi)*cos(-pi/2+alpha);
            ky(ll,kk,jj+199+1)=Kmax*sin(pi/195*(195-weighted_cum(11-jj)))*sin(pi/195*(195-weighted_cum(11-jj))+phi)*cos(-pi/2+alpha);
            kz(ll,kk,jj+199+1)=Kmax*sin(pi/195*(195-weighted_cum(11-jj)))*sin(-pi/2+alpha);
        end
    end
end
kx=permute(kx,[3,2,1]);
ky=permute(ky,[3,2,1]);
kz=permute(kz,[3,2,1]);
Kx_2(:,:,:)=kx([2:106],:,:);%+0.5*kx([1:105],:,:);
Ky_2(:,:,:)=ky([2:106],:,:);%0.5*ky([1:105],:,:);
Kz_2(:,:,:)=kz([2:106],:,:);%0.5*kz([1:105],:,:);

% Kx_2=permute(kx,[2,1]);
Kx_2=Kx_2(:);
Ky_2=Ky_2(:);
Kz_2=Kz_2(:);
Kx_3=pi/max(abs(Kx_2))*Kx_2;
Ky_3=pi/max(abs(Ky_2))*Ky_2;
Kz_3=pi/max(abs(Kz_2))*Kz_2;

mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([Kx_3(:) Ky_3(:) Kz_3(:)]/(4*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(Kx_3,1)*1; % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=1*(ir_mri_density_comp([Kx_3(:) Ky_3(:) Kz_3(:)]/(4*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';

w2 = kdens/max(kdens(:));
Kx_3_new=reshape(Kx_3,105,190,95);
Ky_3_new=reshape(Ky_3,105,190,95);
Kz_3_new=reshape(Kz_3,105,190,95);
Kx_3_new_new=Kx_3_new(1:end,:,:);
Ky_3_new_new=Ky_3_new(1:end,:,:);
Kz_3_new_new=Kz_3_new(1:end,:,:);
traj_uzay(1,:,:)=squeeze(reshape(Kx_3_new_new,105*190,95));
traj_uzay(2,:,:)=squeeze(reshape(Ky_3_new_new,105*190,95));
traj_uzay(3,:,:)=squeeze(reshape(Kz_3_new_new,105*190,95));

traj_rad2 = bart('scale 20.2', traj_uzay);
data_first=inref(3:2:420,:,:);
data_second=inref(4:2:420,:,:);
data_all=(data_first([1:105],:,:)+data_second([1:105],:,:))/2;
clear data_first;
clear data_second;
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ');  
ncoils_compressed = max(find(diag(S)/S(1)>0.05)); 
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
ncoils_compressed

%%
clear inref;
clear data;
[refmrprot, refmdh, inref] =rdMeas_dene('/working/larson2/xin/Revive_03_22_24/meas_MID00129_FID23370_fid_rosette_2mm2x_TE1600.dat');
data_first=inref(3:2:420,:,:);
data_second=inref(4:2:420,:,:);
data_all=(data_first([1:105],:,:)+data_second([1:105],:,:))/2;
clear data_first;
clear data_second;
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ'); 
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
clear x;
for abc=1:ncoils_compressed
      lll=data(:,:,abc);
xcp=(Gm1'*(((lll(:)).*(w2)')));
	xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);

x(:,:,:,abc)=embed(xpcg,mask2(:,:,:));
end

% Coil sensitivity estimation and combination
lowres_img = bart('nufft -i -d30:30:30 -t', traj_rad2, reshape(data,1,190*105,95,ncoils_compressed));
lowres_ksp = bart('fft -u 7', lowres_img);
ksp_zerop = bart('resize -c 0 128 1 128 2 128', lowres_ksp);
sens = bart('ecalib -t 0.0001 -m1', ksp_zerop);
zz1_1600(:,:,:)= sum(conj(sens) .*x ,4) ./ 1;

%%
% for jj=1:128
%     px = angle(squeeze((x(:,:,jj,:))));
% % px = angle(squeeze((fid_data_block1(:,:,2,:,1))));
% nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
% ll=estimate_csm_walsh(nws_water_nuf);
% ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
%     
% nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
% img_subject_5_360_18(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
% end
%%
img_first_echo=zeros(128,128,128,24);
img_first_echo(:,:,:,1)=zz1_20;
img_first_echo(:,:,:,2)=zz1_40;
img_first_echo(:,:,:,3)=zz1_80;
img_first_echo(:,:,:,4)=zz1_160;
img_first_echo(:,:,:,5)=zz1_400;
img_first_echo(:,:,:,6)=zz1_1000;
img_first_echo(:,:,:,7)=zz1_1500;
img_first_echo(:,:,:,8)=zz1_2100;
img_first_echo(:,:,:,9)=zz1_20_12;
img_first_echo(:,:,:,10)=zz1_40_12;
img_first_echo(:,:,:,11)=zz1_80_12;
img_first_echo(:,:,:,12)=zz1_160_12;
img_first_echo(:,:,:,13)=zz1_400_12;
img_first_echo(:,:,:,14)=zz1_1000_12;
img_first_echo(:,:,:,15)=zz1_1500_12;
img_first_echo(:,:,:,16)=zz1_2100_12;
img_first_echo(:,:,:,17)=zz1_20_18;
img_first_echo(:,:,:,18)=zz1_40_18;
img_first_echo(:,:,:,19)=zz1_80_18;
img_first_echo(:,:,:,20)=zz1_160_18;
img_first_echo(:,:,:,21)=zz1_400_18;
img_first_echo(:,:,:,22)=zz1_1000_18;
img_first_echo(:,:,:,23)=zz1_1500_18;
img_first_echo(:,:,:,24)=zz1_2100_18;
%%
clear Kx_2 Ky_2 Kz_2 Kx_3 Ky_3 Kz_3;
% Kx_2(:,:,:)=(kx([105:end-1],:,:));
% Ky_2(:,:,:)=(ky([105:end-1],:,:));
% Kz_2(:,:,:)=(kz([105:end-1],:,:));
Kx_2(:,:,:)=(kx([105:end-1],:,:)+kx([106:end],:,:))/2;
Ky_2(:,:,:)=(ky([105:end-1],:,:)+ky([106:end],:,:))/2;
Kz_2(:,:,:)=(kz([105:end-1],:,:)+kz([106:end],:,:))/2;
Kx_2=Kx_2(:);
Ky_2=Ky_2(:);
Kz_2=Kz_2(:);
Kx_3=(9/9)*pi/max(abs(Kx_2))*Kx_2;
Ky_3=(9/9)*pi/max(abs(Ky_2))*Ky_2;
Kz_3=(9/9)*pi/max(abs(Kz_2))*Kz_2;

Kx_3_new=reshape(Kx_3,105,190,95);
Ky_3_new=reshape(Ky_3,105,190,95);
Kz_3_new=reshape(Kz_3,105,190,95);
Kx_3_new_new=Kx_3_new(1:end,:,:);
Ky_3_new_new=Ky_3_new(1:end,:,:);
Kz_3_new_new=Kz_3_new(1:end,:,:);
traj_uzay(1,:,:)=squeeze(reshape(Kx_3_new_new,105*190,95));
traj_uzay(2,:,:)=squeeze(reshape(Ky_3_new_new,105*190,95));
traj_uzay(3,:,:)=squeeze(reshape(Kz_3_new_new,105*190,95));
traj_rad2 = bart('scale 20.2', traj_uzay);

mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([Kx_3(:) Ky_3(:) Kz_3(:)]/(4*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(Kx_3,1)*1; % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=1*(ir_mri_density_comp([Kx_3(:) Ky_3(:) Kz_3(:)]/(4*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));
data_first=inref(1:2:420,:,:);
data_second=inref(2:2:420,:,:);
data_all=(data_first([106:end],:,:)+data_second([106:end],:,:))/2;
clear data_first;
clear data_second;
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ');  
% Nc = max(find(diag(S)/S(1)>0.01)); 
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
ncoils_compressed
%%
clear inref;
clear data;
[refmrprot, refmdh, inref] =rdMeas_dene('meas_MID00129_FID23370_fid_rosette_2mm2x_TE1600.dat');
data_first=inref(1:2:420,:,:);
data_second=inref(2:2:420,:,:);
data_all=(data_first([106:end],:,:)+data_second([106:end],:,:))/2;
clear data_first;
clear data_second;
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ'); 
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
lowres_img = bart('nufft -i -d30:30:30 -t', traj_rad2, reshape(data,1,190*105,95,ncoils_compressed));
lowres_ksp = bart('fft -u 7', lowres_img);
ksp_zerop = bart('resize -c 0 128 1 128 2 128', lowres_ksp);
sens = bart('ecalib -t 0.0001 -m1', ksp_zerop);
clear x;
for abc=1:ncoils_compressed
      lll=data(:,:,abc);
xcp=(Gm1'*(((lll(:)).*(w2)')));
	xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);

x(:,:,:,abc)=embed(xpcg,mask2(:,:,:));
end 
zz1_second_1600(:,:,:)= sum(conj(sens) .*x ,4) ./ 1;
%%
img_second_echo=zeros(128,128,128,24);
img_second_echo(:,:,:,1)=zz1_2120;
img_second_echo(:,:,:,2)=zz1_2140;
img_second_echo(:,:,:,3)=zz1_2180;
img_second_echo(:,:,:,4)=zz1_2260;
img_second_echo(:,:,:,5)=zz1_2500;
img_second_echo(:,:,:,6)=zz1_3100;
img_second_echo(:,:,:,7)=zz1_3600;
img_second_echo(:,:,:,8)=zz1_4200;
img_second_echo(:,:,:,9)=zz1_2120_12;
img_second_echo(:,:,:,10)=zz1_2140_12;
img_second_echo(:,:,:,11)=zz1_2180_12;
img_second_echo(:,:,:,12)=zz1_2260_12;
img_second_echo(:,:,:,13)=zz1_2500_12;
img_second_echo(:,:,:,14)=zz1_3100_12;
img_second_echo(:,:,:,15)=zz1_3600_12;
img_second_echo(:,:,:,16)=zz1_4200_12;
img_second_echo(:,:,:,17)=zz1_2120_18;
img_second_echo(:,:,:,18)=zz1_2140_18;
img_second_echo(:,:,:,19)=zz1_2180_18;
img_second_echo(:,:,:,20)=zz1_2260_18;
img_second_echo(:,:,:,21)=zz1_2500_18;
img_second_echo(:,:,:,22)=zz1_3100_18;
img_second_echo(:,:,:,23)=zz1_3600_18;
img_second_echo(:,:,:,24)=zz1_4200_18;