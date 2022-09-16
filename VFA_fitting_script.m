%% test processing script to combine multiple flip angle datasets

close all;

%% Load data
% First load the Pfiles for the 18, 12, 6 and AFI scans

addpath(genpath('/home/sf719662/Documents/UTEMRI_Brain'));
addpath(genpath('/data/larson2/brain_uT2/orchestra-sdk-1.7-1.matlab'));
% addpath(genpath('/data/larson2/brain_uT2/2019-08-01_3T_8TE'));
datapath = '/data/larson2/brain_uT2/';
suffix = '-csreconallec_l2_r0p01.mat';
scan_date = '2020-01-28';
dateFolder = '_3T_8TE/';
Pfile_18deg = 'P50176.7_01281445';
Pfile_12deg = 'P51712.7_01281458';
Pfile_06deg = 'P53248.7_01281512';
Pfile_AFI = 'P54272.7_01281523';




% In the case of unequal R1 and R2 values, we can read in auto prescan values to correct scaling
pfile_18 = GERecon('Pfile.Load', Pfile_18deg);
header_18 = GERecon('Pfile.Header', pfile_18);

pfile_12 = GERecon('Pfile.Load', Pfile_12deg);
header_12 = GERecon('Pfile.Header', pfile_12);

pfile_06 = GERecon('Pfile.Load', Pfile_06deg);
header_06 = GERecon('Pfile.Header', pfile_06);


% scaling correction if r1 and/or r2 do not match between scans normalized to 18 degree flip scan 
% unequal values for TG cannot be corrected for analytically
auto_prescan_values = struct('aps_r1', {}, 'aps_r2', {}, 'aps_tg', {});
auto_prescan_values(1).aps_r1 = [header_18.PrescanHeader.aps_r1, header_12.PrescanHeader.aps_r1, header_06.PrescanHeader.aps_r1];
auto_prescan_values(1).aps_r2 = [header_18.PrescanHeader.aps_r2, header_12.PrescanHeader.aps_r2, header_06.PrescanHeader.aps_r2];
auto_prescan_values(1).aps_tg = [header_18.PrescanHeader.aps_tg, header_12.PrescanHeader.aps_tg, header_06.PrescanHeader.aps_tg];


% 18-deg
load([datapath scan_date dateFolder Pfile_18deg suffix]);

imall_all{3} = imall ; 
imallplus_all{3} = imallplus;
TEin{3} = TE;
flips(3) = 18*pi/180;

% 12-deg
load([datapath scan_date dateFolder Pfile_12deg suffix]);

imall_all{2} = imall .* (10 ^ (3 * (header_18.PrescanHeader.aps_r1 - header_12.PrescanHeader.aps_r1) / 20)) .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_12.PrescanHeader.aps_r2)); 
imallplus_all{2} = imallplus .* (10 ^ (3 * (header_18.PrescanHeader.aps_r1 - header_12.PrescanHeader.aps_r1) / 20)) .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_12.PrescanHeader.aps_r2)); 
flips(2) = 12*pi/180;
TEin{2} = TE;

% 6-deg
load([datapath scan_date dateFolder Pfile_06deg suffix]);

imall_all{1} = imall .* 10 ^(3 * (header_18.PrescanHeader.aps_r1 - header_06.PrescanHeader.aps_r1) /20) .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_06.PrescanHeader.aps_r2)); 
imallplus_all{1} = imallplus .* 10 ^(3 * (header_18.PrescanHeader.aps_r1 - header_06.PrescanHeader.aps_r1) /20) .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_06.PrescanHeader.aps_r2)); 
flips(1) = 6*pi/180;
TEin{1} = TE;

for j = 1:3
    TEin{j} = TEin{j} * 1e-3;
end

TR = 9.1 * 1e-3; % s

% Visualize flip angle volumes if necessary
figure(1)
subplot(131)
imagesc(flip(abs(imallplus_all{3}(:,:,46)))); colormap gray; axis off; colorbar
subplot(132)
imagesc(flip(abs(imallplus_all{2}(:,:,46))));colormap gray; axis off; colorbar
subplot(133)
imagesc(flip(abs(imallplus_all{1}(:,:,46))));colormap gray; axis off; colorbar


%% B1 correction

% Now we'll load in the AFI 

AFI_nominal_flip = 45 * pi/180; % UTE AFIcones

load([Pfile_AFI '_UTE_AFI.mat']); % load in UTE AFI volumes and extrac the 2 volumes
AFI_map_1 = squeeze(recon_grid(:,:,1,:));
AFI_map_2 = squeeze(recon_grid(:,:,2,:));

% Generate B1+ map analytically as described in Yarnykh 2006.
S1 = abs(AFI_map_1); S2 = abs(AFI_map_2); % magnitude images
TR1 = 7; TR2 = 35; % ms
n = TR2 / TR1; r = S2 ./ S1;
FA_map = abs(acos((n * r - 1) ./ (n - r))); 

FA_map_resize = FA_map(:, 6:91, :); % make sure dimensions match UTE image dimensions

% Now we will match the orientation of the B1+ map to match the UTE scans
FA_map_ute_orientation = flip(flip(FA_map_resize, 3), 1); % 3D UTE AFI cones

%% Visualize B1+ map orientation

% We can visualize the B1+ map if desired

figure(1) % UTE scan orientation
subplot(131)
imagesc(abs(flip(imall_all{1}(:,:,46)))); colormap gray; 
subplot(132)
imagesc(imrotate(abs(squeeze(imall_all{1}(:,50,:,1))), 270)); colormap gray
subplot(133)
imagesc(imrotate(abs(squeeze(imall_all{1}(50,:,:,1))), 270)); colormap gray

figure(2) % B1+ map orientation
%subplot(131)
imagesc(flip(FA_map_ute_orientation(:,:,46)), [0 1]); colormap gray; axis off;
%subplot(132)
imagesc(flip(imrotate(squeeze(FA_map_ute_orientation(:,50,:,1)), 270), 2), [0 1]); colormap gray; axis off
%subplot(133)
imagesc(imrotate(squeeze(FA_map_ute_orientation(50,:,:,1)), 270), [0 1]); colormap gray; axis off; colorbar




%% brainmask

% Now we will generate a brainmask that isolates the cerebrum and removes
% extraneous anatomical regions that we do not wish to fit to. This will speed up computation time.

for n = 1:3 %3
    im_diff{n} = (abs(imallplus_all{n}(:,:,:,1)) - abs(imallplus_all{n}(:,:,:,8)));

    brainmask{n} = (im_diff{n} < 0.034) & abs(imallplus_all{n}(:,:,:,1)) > mean(mean(mean(abs(imallplus_all{n}(:,:,:,1)))))*.5;

    se = strel('disk', 5);
    brainmask{n} = imopen(brainmask{n}, se);

    brainmask{n}(54:end,:,60:end) = 0; 
    brainmask{n}(73:end,:,55:end) = 0; 
    brainmask{n}(95:end,:,44:end) = 0; 
    brainmask{n}(:,:,73:end) = 0;  
    brainmask{n}(:,:,1:10) = 0; 
    % 
    im_masked{n} = imallplus_all{n}(:,:,:,1) .* brainmask{n};
end

% We can visualize the brainmask in 3 orientations if desired
figure
subplot(131)
imagesc(squeeze(abs(flip(im_masked{3}(:,:,46,1))))); colormap gray
subplot(132)
imagesc(imrotate(squeeze(abs(flip(im_masked{3}(:,50,:,1)))), 270));
subplot(133)
imagesc(imrotate(squeeze(abs(flip(im_masked{3}(50,:,:,1)))), 270));

%% setup data
% Now we will prepare the data to be fed into the fitting model

B0 = 3;
phi_RF = .0562;
fit_thresh = [.011 .025];


n=2; % 2
    if exist('brainmask')
        im_uteplus_sample = im_masked{n}(:,:,:,1);            % brainmask passed into fitting
    else
        im_uteplus_sample = imallplus_all{n}(:,:,:,1);
    end
    
fit_mask= repmat(abs(im_uteplus_sample) > fit_thresh(1), [1,1,1,size(imallplus_all{n},4)]);

for n=1:3 %3
    imallplus_all{n} = imallplus_all{n} .* fit_mask;
    im_uteplus{n} = imallplus_all{n}(:,:,:,1) .* fit_mask(:,:,:,1);
end


% get size and number of echoes
imsize = size(im_uteplus_sample);

Nechoes = size(imallplus_all{1}, 4);

% only fit non-zero (can be used for masking) values
I = find(im_uteplus{1});

NI = length(I);
        
% get max voxel value
Snormall = max([im_uteplus{1}(:);im_uteplus{2}(:);im_uteplus{3}(:)]);

imfit = cell(1,3); 
imfitt = cell(1,3); 

for n=1:3 % 3
       imfit{n} = reshape(imallplus_all{n}, [prod(imsize) Nechoes]) / Snormall; 
       imfitt{n} = imfit{n}(I,:);
end           


% throw out any nan values and generate positive non-zero AFI volume

if exist('FA_map_ute_orientation')
    if size(find(FA_map_ute_orientation(I)), 1) < size(I, 1)
        x = FA_map_ute_orientation; 
        x(x <= 0) = AFI_nominal_flip;
        x(isnan(x)) = AFI_nominal_flip;
        nzFA_map = x;        
    end
    nzFA_map = FA_map_ute_orientation;
    imagesc(flip(nzFA_map(:,:,46))); colorbar % visualize non-zero FA map if desired
end


%% Fitting_multiscan_test
% Now we can perform multi-component fitting for each flip angle

clear('fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3');
% make sure subet of voxels for fitting is all the same 
parfor Ix = 1:NI
    
    Sin = cell(1,3);
    for n = 1:3 
      Sin{n} = imfitt{n}(Ix,:); 
    end
    
    [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
        utebrain_multiscan_fitting_function(TEin, Sin, B0, phi_RF, 0);
end


savefname = [scan_date '_multiscan-fitting_plus'];
writedata_flag = 1;

if writedata_flag
    save(savefname,'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
end   


fit_filename = [datapath scan_date dateFolder filesep savefname];
fit_plus_filename = [datapath scan_date dateFolder filesep savefname]; 
AIC_thresh = 20;
T1_flag = 0;

% Once fitting is complete we can generate fit maps
fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh);
save('fit_maps_multiscan', 'fit_maps')


%% fitting_T1_test
% Here we can perform multi-component fitting with T1 SPGR term

clear('fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3');

% make sure subet of voxels for fitting is all the same 
parfor Ix = 1:NI
    
    Sin = cell(1,3);
    for n = 1:3
      Sin{n} = imfitt{n}(Ix,:); 
    end
    
    % Calculate actual flip angle per voxel
    B1_scale_map  = nzFA_map(I(Ix)) ./ AFI_nominal_flip;
    aflips = flips .* B1_scale_map;
    
    [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
        utebrain_t1_fitting_function(TEin, Sin, aflips, TR, B0, phi_RF, 1);
end

if exist('nzFA_map')
    savefname = [scan_date '_T1-fitting_lb01_plus_aps_scale_TE_test_B1cor_NEW'];
else
    savefname = [scan_date '_T1-fitting_lb01_plus_aps_scale_TE_test_NEW'];
end
    
writedata_flag = 1;

if writedata_flag
    save(savefname,'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
end

% generate fit maps
fit_filename = [datapath scan_date dateFolder filesep savefname];
fit_plus_filename = [datapath scan_date dateFolder filesep savefname]; 
AIC_thresh = 20;
T1_flag = 1;

fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh);
save('fit_maps_aps_scaled_TE_test_B1cor_NEW', 'fit_maps')

return

    