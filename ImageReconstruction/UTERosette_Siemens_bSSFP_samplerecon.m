data_file = '/working/larson3/Xin/VFA_MS_1/meas_MID03285_FID24049_fid_rosette_2mm_bSSFP_41k.dat';
%data_file = '/working/larson3/Xin/VFA_bssfp_5/meas_MID01732_FID26323_fid_rosette_2mm_bSSFP_41k.dat';
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

%% scan and recon parameters
gamma=42.576; %kHZ/mT
TR=2500; %us
sampling_interval=10; %us
FOV=0.24; %m
T=1.9; %ms
matrix_size=240; % for trajectory

% for nufft
nx=256;
ny=256;
nz=256;

ncoils  =   size(inref,3);
nt=1;
npetals = size(inref,2);
nsamples = size(inref,1);


%% Data modulations

% phase correction for chopping(?)
inref(:,2:2:end,:) = -inref(:,2:2:end,:);

f_offset = 300; %Hz

t=[0:nsamples-1] * sampling_interval/2;  % us
f_modulation = exp(1j*2*pi*f_offset*t'/1e6);
inref = inref .* repmat(f_modulation, [1, npetals, ncoils]);

%% k-space trajectory

Kmax=matrix_size/FOV/2;
K_interval=Kmax*2/matrix_size;
n=1:1:60;
kr=n*K_interval;

Gr=kr*2*pi/(gamma*T);
dead_time=0.2;%ms
Gx_initial=2*kr/gamma/dead_time;
Gy_initial=-Gr;
time_pre=0.2;%ms
weighted=[0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95];
weighted_cum=cumsum(weighted);
for ll=1:190
    alpha=pi/189*(ll-1);
    for kk=1:378
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
        Sx(ll,kk,1)=0;
        Sy(ll,kk,1)=0;
        Sz(ll,kk,1)=0;
        phi=2*pi/378*(kk-1);
        phi_next=2*pi/190*kk;
        for jj=1:20
            kx(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*cos(pi/195*weighted_cum(jj)+phi)*cos(-pi/2+alpha);
            ky(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*sin(pi/195*weighted_cum(jj)+phi)*cos(-pi/2+alpha);
            kz(ll,kk,jj+2)=Kmax*sin(pi/195*weighted_cum(jj))*sin(-pi/2+alpha);
            %     Gx(ll,kk,jj+2)=(kx(ll,kk,jj+2)-kx(ll,kk,jj+1))/gamma/0.01;
            %     Gy(ll,kk,jj+2)=(ky(ll,kk,jj+2)-ky(ll,kk,jj+1))/gamma/0.01;
            %     Gz(ll,kk,jj+2)=(kz(ll,kk,jj+2)-kz(ll,kk,jj+1))/gamma/0.01;
            %     Sx(ll,kk,jj+1)=(Gx(ll,kk,jj+2)-Gx(ll,kk,jj+1))/0.01;
            %     Sy(ll,kk,jj+1)=(Gy(ll,kk,jj+2)-Gy(ll,kk,jj+1))/0.01;
            %     Sz(ll,kk,jj+1)=(Gz(ll,kk,jj+2)-Gz(ll,kk,jj+1))/0.01;
        end
        for ii=11:184
            kx(ll,kk,ii+11+1)=Kmax*sin(pi/195*ii)*cos(pi/195*ii+phi)*cos(-pi/2+alpha);
            ky(ll,kk,ii+11+1)=Kmax*sin(pi/195*ii)*sin(pi/195*ii+phi)*cos(-pi/2+alpha);
            kz(ll,kk,ii+11+1)=Kmax*sin(pi/195*ii)*sin(-pi/2+alpha);
            %     Gx(ll,kk,ii+11+1)=(kx(ll,kk,ii+11+1)-kx(ll,kk,ii+11))/gamma/0.01;
            %     Gy(ll,kk,ii+11+1)=(ky(ll,kk,ii+11+1)-ky(ll,kk,ii+11))/gamma/0.01;
            %     Gz(ll,kk,ii+11+1)=(kz(ll,kk,ii+11+1)-kz(ll,kk,ii+11))/gamma/0.01;
            %     Sx(ll,kk,ii+11)=(Gx(ll,kk,ii+11+1)-Gx(ll,kk,ii+11))/0.01;
            %     Sy(ll,kk,ii+11)=(Gy(ll,kk,ii+11+1)-Gy(ll,kk,ii+11))/0.01;
            %     Sz(ll,kk,ii+11)=(Gz(ll,kk,ii+11+1)-Gz(ll,kk,ii+11))/0.01;
        end
        for jj=1:20
            kx(ll,kk,jj+195+1)=Kmax*sin(pi/195*(195-weighted_cum(21-jj)))*cos(pi/195*(195-weighted_cum(21-jj))+phi)*cos(-pi/2+alpha);
            ky(ll,kk,jj+195+1)=Kmax*sin(pi/195*(195-weighted_cum(21-jj)))*sin(pi/195*(195-weighted_cum(21-jj))+phi)*cos(-pi/2+alpha);
            kz(ll,kk,jj+195+1)=Kmax*sin(pi/195*(195-weighted_cum(21-jj)))*sin(-pi/2+alpha);
            %     Gx(ll,kk,jj+195+1)=(kx(ll,kk,jj+195+1)-kx(ll,kk,jj+195))/gamma/0.01;
            %     Gy(ll,kk,jj+195+1)=(ky(ll,kk,jj+195+1)-ky(ll,kk,jj+195))/gamma/0.01;
            %     Gz(ll,kk,jj+195+1)=(kz(ll,kk,jj+195+1)-kz(ll,kk,jj+195))/gamma/0.01;
            %     Sx(ll,kk,jj+195)=(Gx(ll,kk,jj+195+1)-Gx(ll,kk,jj+195))/0.01;
            %     Sy(ll,kk,jj+195)=(Gy(ll,kk,jj+195+1)-Gy(ll,kk,jj+195))/0.01;
            %     Sz(ll,kk,jj+195)=(Gz(ll,kk,jj+195+1)-Gz(ll,kk,jj+195))/0.01;
        end
    end
end
kx=permute(kx,[3,2,1]);
ky=permute(ky,[3,2,1]);
kz=permute(kz,[3,2,1]);
% Kx_2(:,:,:)=0.5*kx([2:106],:,:)+0.5*kx([1:105],:,:);
% Ky_2(:,:,:)=0.5*ky([2:106],:,:)+0.5*ky([1:105],:,:);
% Kz_2(:,:,:)=0.5*kz([2:106],:,:)+0.5*kz([1:105],:,:);
%%
Kx_new=zeros(432,378,190);
Ky_new=zeros(432,378,190);
Kz_new=zeros(432,378,190);
for ii=1:431
    if rem(ii,2)==1
        Kx_new(ii,:,:)=kx(fix(ii/2)+1,:,:);
        Ky_new(ii,:,:)=ky(fix(ii/2)+1,:,:);
        Kz_new(ii,:,:)=kz(fix(ii/2)+1,:,:);
    else
        Kx_new(ii,:,:)=(kx(fix(ii/2),:,:)+kx(fix(ii/2)+1,:,:))/2;
        Ky_new(ii,:,:)=(ky(fix(ii/2),:,:)+ky(fix(ii/2)+1,:,:))/2;
        Kz_new(ii,:,:)=(kz(fix(ii/2),:,:)+kz(fix(ii/2)+1,:,:))/2;
    end
end
Kx_new(432,:,:)=kx(216,:,:);
Ky_new(432,:,:)=ky(216,:,:);
Kz_new(432,:,:)=kz(216,:,:);
% %%
% Kx_new_2=zeros(216,378,190);
% Ky_new_2=zeros(216,378,190);
% Kz_new_2=zeros(216,378,190);
% for ii=1:216
%     for jj=1:378
%         if rem(jj,2)==1
%             Kx_new_2(ii,jj,:)=Kx_new(2*(ii-1)+1,jj,:);
%             Ky_new_2(ii,jj,:)=Ky_new(2*(ii-1)+1,jj,:);
%             Kz_new_2(ii,jj,:)=Kz_new(2*(ii-1)+1,jj,:);
%         else
%             Kx_new_2(ii,jj,:)=Kx_new(2*(ii),jj,:);
%             Ky_new_2(ii,jj,:)=Ky_new(2*(ii),jj,:);
%             Kz_new_2(ii,jj,:)=Kz_new(2*(ii),jj,:);
%         end
%     end
% end

%%
% Kx_2=permute(kx,[2,1]);
Kx_2(:,:,:)=Kx_new([2:217],:,:);
Ky_2(:,:,:)=Ky_new([2:217],:,:);
Kz_2(:,:,:)=Kz_new([2:217],:,:);
% Kx_2(:,:,:)=Kx_new_2([2:109],:,:);
% Ky_2(:,:,:)=Ky_new_2([2:109],:,:);
% Kz_2(:,:,:)=Kz_new_2([2:109],:,:);
% time=40:5:215*5+40;
% decay=exp(-time/500);
% % fun=@(x) exp(-x/500);
% decay=decay';
% decay=repmat(decay,1,378,190);
Kx_2=Kx_2(:);
Ky_2=Ky_2(:);
Kz_2=Kz_2(:);
Kx_3=pi/max(abs(Kx_2))*Kx_2;
Ky_3=pi/max(abs(Ky_2))*Ky_2;
Kz_3=pi/max(abs(Kz_2))*Kz_2;

%% Reconstruction

mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([Kx_3(:) Ky_3(:) Kz_3(:)]/(2*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(Kx_3,1)*1; % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=1*(ir_mri_density_comp([Kx_3(:) Ky_3(:) Kz_3(:)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));
% w2 = kdens/sum(kdens(:));
% kspace_simu_data=ones(15513120,1);
% xcp=(Gm1'*(((kspace_simu_data).*(w2)'.*decay(:))));
% 	xpcg = qpwls_pcg1(xcp, Gm1, 1, kspace_simu_data.*decay(:), R.C, 'niter', 1);
% 
% x(:,:,:)=embed(xpcg,mask2(:,:,:));
% fwhm = qpwls_psf_fwhm(x, ndim, arg, units);
% w2 = decay(:)'.*kdens/max(kdens(:));
% zzzz=qpwls_psf(Gm1,R,1,mask2,w2.','fwhmtype','profile');
% zzzz=qpwls_psf(Gm1,R,1,mask2,decay(:),'fwhmtype','profile');
Kx_3_new=reshape(Kx_3,216,378,190);
Ky_3_new=reshape(Ky_3,216,378,190);
Kz_3_new=reshape(Kz_3,216,378,190);
Kx_3_new_new=Kx_3_new(1:end,:,:);
Ky_3_new_new=Ky_3_new(1:end,:,:);
Kz_3_new_new=Kz_3_new(1:end,:,:);
traj_uzay(1,:,:)=squeeze(reshape(Kx_3_new_new,216*378,190));
traj_uzay(2,:,:)=squeeze(reshape(Ky_3_new_new,216*378,190));
traj_uzay(3,:,:)=squeeze(reshape(Kz_3_new_new,216*378,190));
% traj_rad2 = bart('scale 20.2', traj_uzay);
% data_first=inref(3:2:420,:,:);
% data_second=inref(4:2:420,:,:);
% data_all=inref_undersample(2:109,:,:);
data_all=inref(3:218,:,:);
clear data_first;
clear data_second;
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ');  
ncoils_compressed = max(find(diag(S)/S(1)>0.01)); %0.01
%%
% clear inref;
% clear data;
% clear D;
% clear V;
% clear img_subject_5;
% [refmrprot, refmdh, inref] =rdMeas_dene('/working/larson3/Xin/VFA_bssfp_5_repeat/meas_MID02768_FID59328_fid_rosette_2mm_bSSFP_41k.dat');
% for ii=2:2:71820
%     for jj=1:58
%         inref(:,ii,jj)=inref(:,ii,jj).*exp(1j*pi);
%     end
% end
% t=0:5:431*5;
% t=t/1000;
% for ii=1:1:71820
%     for jj=1:58
%         inref(:,ii,jj)=inref(:,ii,jj).*exp(1j*2*pi*300*t'/1000);
%     end
% end
% data_all=inref(3:218,:,:);
% % data_all=inref(2:217,:,:);
% clear data_first;
% clear data_second;
% data=squeeze(data_all); 
% clear data_all;
% D = reshape(data,size(data,1)*size(data,2),nc);
% [U,S,V] = svd(D,'econ');
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
% lowres_img = bart('nufft -i -d30:30:30 -t', traj_rad2, reshape(data,1,216*378,190,ncoils_compressed));
% lowres_ksp = bart('fft -u 7', lowres_img);
% ksp_zerop = bart('resize -c 0 256 1 256 2 256', lowres_ksp);
% sens = bart('ecalib -t 0.0001 -m1', ksp_zerop);
% zz1_2(:,:,:)= sum(conj(sens) .*x ,4) ./ 1;
for jj=1:256
    px = angle(squeeze((x(:,:,jj,:))));
% px = angle(squeeze((fid_data_block1(:,:,2,:,1))));
nws_water_nuf=squeeze((x(:,:,jj,:))).* exp( -1i * px );
ll=estimate_csm_walsh(nws_water_nuf);
ll_sq = sum(ll .* conj(ll),3); ll(ll < eps) = 1;
    
nws_water_nuf=squeeze((x(:,:,jj,:))/1).* exp( -1i * px );
img_recon(:,:,jj)= sum(conj(ll) .*nws_water_nuf ,3) ./ ll_sq;
end

%%  second TE
clear Kx_2 Ky_2 Kz_2 Kx_3 Ky_3 Kz_3 data;
% Kx_2(:,:,:)=(kx([105:end-1],:,:));
% Ky_2(:,:,:)=(ky([105:end-1],:,:));
% Kz_2(:,:,:)=(kz([105:end-1],:,:));
Kx_2(:,:,:)=0.35*Kx_new([215:430],:,:)+0.65*Kx_new([216:431],:,:);
Ky_2(:,:,:)=0.35*Ky_new([215:430],:,:)+0.65*Ky_new([216:431],:,:);
Kz_2(:,:,:)=0.35*Kz_new([215:430],:,:)+0.65*Kz_new([216:431],:,:);
% Kx_2(:,:,:)=0.5*Kx_new([215:430],:,:)+0.5*Kx_new([216:431],:,:);
% Ky_2(:,:,:)=0.5*Ky_new([215:430],:,:)+0.5*Ky_new([216:431],:,:);
% Kz_2(:,:,:)=0.5*Kz_new([215:430],:,:)+0.5*Kz_new([216:431],:,:);
Kx_2=Kx_2(:);
Ky_2=Ky_2(:);
Kz_2=Kz_2(:);
Kx_3=(9/9)*pi/max(abs(Kx_2))*Kx_2;
Ky_3=(9/9)*pi/max(abs(Ky_2))*Ky_2;
Kz_3=(9/9)*pi/max(abs(Kz_2))*Kz_2;

Kx_3_new=reshape(Kx_3,216,378,190);
Ky_3_new=reshape(Ky_3,216,378,190);
Kz_3_new=reshape(Kz_3,216,378,190);
Kx_3_new_new=Kx_3_new(1:end,:,:);
Ky_3_new_new=Ky_3_new(1:end,:,:);
Kz_3_new_new=Kz_3_new(1:end,:,:);
traj_uzay(1,:,:)=squeeze(reshape(Kx_3_new_new,216*378,190));
traj_uzay(2,:,:)=squeeze(reshape(Ky_3_new_new,216*378,190));
traj_uzay(3,:,:)=squeeze(reshape(Kz_3_new_new,216*378,190));
% traj_rad2 = bart('scale 20.2', traj_uzay);
mask2 = true([nx ny nz]);
nufft2 = {[nx ny nz], [3 3 3], 2*[nx ny nz], [nx/2 ny/2 nz/2 ], 'table', 2^10, 'minmax:kb'};
f.basis = {'rect'};
Gm1 = Gmri([Kx_3(:) Ky_3(:) Kz_3(:)]/(2*pi), mask2, 'fov', 256, 'basis', f.basis, 'nufft', nufft2);
beta = 2^-21 * size(Kx_3,1)*1; % good for quadratic 'rect'
R = Reg1(mask2, 'beta', beta);
kdens=1*(ir_mri_density_comp([Kx_3(:) Ky_3(:) Kz_3(:)]/(2*pi),'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',256}))';
w2 = kdens/max(kdens(:));
data_all=inref([217:432],:,:);
% data_all=inref([216:431],:,:);
data=squeeze(data_all); 
clear data_all;
D = reshape(data,size(data,1)*size(data,2),ncoils);
[U,S,V] = svd(D,'econ');  
% ncoils_compressed = max(find(diag(S)/S(1)>0.01)); 
data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
%%
% clear inref;
% clear data;
% clear D;
% clear V;
% clear img_subject;
% clear     second_5;
% [refmrprot, refmdh, inref] =rdMeas_dene('/working/larson3/Xin/VFA_bssfp_5_repeat/meas_MID02768_FID59328_fid_rosette_2mm_bSSFP_41k.dat');
% for ii=2:2:71820
%     for jj=1:58
%         inref(:,ii,jj)=inref(:,ii,jj).*exp(1j*pi);
%     end
% end
% t=0:5:431*5;
% t=t/1000;
% for ii=1:1:71820
%     for jj=1:58
%         inref(:,ii,jj)=inref(:,ii,jj).*exp(1j*2*pi*300*t'/1000);
%     end
% end
% data_all=inref([217:432],:,:);
% % data_all=inref(2:217,:,:);
% clear data_first;
% clear data_second;
% data=squeeze(data_all); 
% clear data_all;
% D = reshape(data,size(data,1)*size(data,2),nc);
% [U,S,V] = svd(D,'econ');
% data = reshape(D*V(:,1:ncoils_compressed),size(data,1),size(data,2),ncoils_compressed);
ncoils_compressed
clear x;
for abc=1:ncoils_compressed
      lll=data(:,:,abc);
xcp=(Gm1'*(((lll(:)).*(w2)')));
	xpcg = qpwls_pcg1(xcp, Gm1, 1, lll(:), R.C, 'niter', 1);

x(:,:,:,abc)=embed(xpcg,mask2(:,:,:));
abc
end
% lowres_img = bart('nufft -i -d30:30:30 -t', traj_rad2, reshape(data,1,216*378,190,ncoils_compressed));
% lowres_ksp = bart('fft -u 7', lowres_img);
% ksp_zerop = bart('resize -c 0 256 1 256 2 256', lowres_ksp);
% sens = bart('ecalib -t 0.0001 -m1', ksp_zerop);
% zz1_2_second(:,:,:)= sum(conj(sens) .*x ,4) ./ 1;
for jj=1:256 
    px = angle(squeeze((x(:,:,jj,:))));
% px = angle(squeeze((fid_data_block1(:,:,2,:,1))));
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