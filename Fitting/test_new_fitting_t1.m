%% test processing script to combine multiple flip angle datasets

close all;

%% load data

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


% read in auto prescan values to correct scaling
pfile_18 = GERecon('Pfile.Load', Pfile_18deg);
header_18 = GERecon('Pfile.Header', pfile_18);

pfile_12 = GERecon('Pfile.Load', Pfile_12deg);
header_12 = GERecon('Pfile.Header', pfile_12);

pfile_06 = GERecon('Pfile.Load', Pfile_06deg);
header_06 = GERecon('Pfile.Header', pfile_06);

auto_prescan_values = struct('aps_r1', {}, 'aps_r2', {}, 'aps_tg', {});
auto_prescan_values(1).aps_r1 = [header_18.PrescanHeader.aps_r1, header_12.PrescanHeader.aps_r1, header_06.PrescanHeader.aps_r1];
auto_prescan_values(1).aps_r2 = [header_18.PrescanHeader.aps_r2, header_12.PrescanHeader.aps_r2, header_06.PrescanHeader.aps_r2];
auto_prescan_values(1).aps_tg = [header_18.PrescanHeader.aps_tg, header_12.PrescanHeader.aps_tg, header_06.PrescanHeader.aps_tg];

% 18-deg
load([datapath scan_date dateFolder Pfile_18deg suffix]);

imall_all{3} = imall; 
imallplus_all{3} = imallplus;
TEin{3} = TE;
flips(3) = 18*pi/180;

% 12-deg
load([datapath scan_date dateFolder Pfile_12deg suffix]);

imall_all{2} = imall;% .* (auto_prescan_values(1).aps_r1(1) / auto_prescan_values(1).aps_r1(2)) .* (auto_prescan_values(1).aps_r2(1) / auto_prescan_values(1).aps_r2(2));  
imallplus_all{2} = imallplus;% .* (auto_prescan_values(1).aps_r1(1) / auto_prescan_values(1).aps_r1(2)) .* (auto_prescan_values(1).aps_r2(1) / auto_prescan_values(1).aps_r2(2)); % testing scale factor
flips(2) = 12*pi/180;
TEin{2} = TE;

% 6-deg
load([datapath scan_date dateFolder Pfile_06deg suffix]);

imall_all{1} = imall;% .* (auto_prescan_values(1).aps_r1(1) / auto_prescan_values(1).aps_r1(3)) .* (auto_prescan_values(1).aps_r2(1) / auto_prescan_values(1).aps_r2(3)); 
imallplus_all{1} = imallplus;% .* (auto_prescan_values(1).aps_r1(1) / auto_prescan_values(1).aps_r1(3)) .* (auto_prescan_values(1).aps_r2(1) / auto_prescan_values(1).aps_r2(3));
flips(1) = 6*pi/180;
TEin{1} = TE;

for j = 1:3
    TEin{j} = TEin{j} * 1e-3;
end

TR = 9.1 * 1e-3; % s

% visualize flip angle volumes and check scaling
figure
subplot(131)
imagesc(flip(abs(imallplus_all{3}(:,:,46)))); colormap gray; colorbar
subplot(132)
imagesc(flip(abs(imallplus_all{2}(:,:,46))));colormap gray; colorbar
subplot(133)
imagesc(flip(abs(imallplus_all{1}(:,:,46))));colormap gray; colorbar

%% B1 correction

% AFI 
%AFI_nominal_flip = 44 * pi/180; % UTE AFI radial
AFI_nominal_flip = 45 * pi/180; % UTE AFIcones


% S1
info_1 = dir('/data/larson2/brain_uT2/2019-11-21_3T_8TE/E8348/6/*DCM'); % MAKE SURE dicoms 1-9 are labeled with 01-09 instead
%info_1 = dir('/data/larson2/brain_uT2/2019-10-17_3T_8TE/Dicom_S8/*dcm');
%info_1 = dir('/data/larson2/brain_uT2/2019-11-21_3T_8TE/E8348/6/*DCM');
% imall_size = size(imall);
% xLength = imall_size(1);
% yLength = imall_size(2);

AFI_map = zeros(96, 96, length(info_1));
for n = 1:length(info_1)
    AFI_map(:,:,n) = dicomread(info_1(n).name);
end

% S2
info_2 = dir('/data/larson2/brain_uT2/2019-11-21_3T_8TE/E8348/6_2/*DCM'); % MAKE SURE dicoms 1-9 are labeled with 01-09 instead
%info_2 = dir('/data/larson2/brain_uT2/2019-10-17_3T_8TE/Dicom_S8_b/*dcm');
%info_2 = dir('/data/larson2/brain_uT2/2019-11-21_3T_8TE/E8348/6_2/*DCM');
% imall_size = size(imall);
% xLength = imall_size(1);
% yLength = imall_size(2);
AFI_map_2nd = zeros(96, 96, length(info_2));
for n = 1:length(info_2)
    AFI_map_2nd(:,:,n) = dicomread(info_2(n).name);
end

% generate flip angle map - AFI
S1 = abs(AFI_map); S2 = abs(AFI_map_2nd); % magnitude images
TR1 = 7; TR2 = 35; % ms
n = TR2 / TR1; r = S2 ./ S1;
FA_map = abs(acos((n * r - 1) ./ (n - r))); 

FA_map_resize = FA_map(:, 6:91, :); % make sure dimensions match UTE image dimensions
% FA_map_deg = FA_map_resize * 180/pi; % convert to degrees if necessary
% FA_map_scaled = FA_map_resize ./ AFI_nominal_flip;
imagesc(FA_map_resize(:,:,46)); colormap gray; colorbar

% hard code correction test
FA_map_resize_cor = FA_map_resize + 0.35;


%% AFI B1 map orientation
%FA_map_ute_orientation = flip(flip(imrotate(FA_map_resize_cor, 270), 2), 3); % 3D UTE AFI radial
% FA_map_ute_orientation = flip(FA_map_resize, 3);
FA_map_ute_orientation = flip(flip(FA_map_resize_cor, 3), 1); % 3D UTE AFI cones


% check orientation
figure(1)
subplot(231)
imagesc(abs(flip(imall_all{1}(:,:,46)))); colormap gray; 
subplot(232)
imagesc(imrotate(abs(squeeze(imall_all{1}(:,50,:,1))), 270)); colormap gray
subplot(233)
imagesc(imrotate(abs(squeeze(imall_all{1}(50,:,:,1))), 270)); colormap gray

subplot(234)
imagesc(flip(FA_map_ute_orientation(:,:,46)), [0 1]); colormap gray;
subplot(235)
imagesc(imrotate(squeeze(FA_map_ute_orientation(:,50,:,1)), 270), [0 1]); colormap gray
subplot(236)
imagesc(imrotate(squeeze(FA_map_ute_orientation(50,:,:,1)), 270), [0 1]); colormap gray




%% brainmask

% im_diff = cell(1, 3);
% 
% if length(im_diff) == 1
%     imall_all{1} = imall;
%     imallplus_all{1} = imallplus;
% end

for n = 1:3 %3
    im_diff{n} = (abs(imallplus_all{n}(:,:,:,1)) - abs(imallplus_all{n}(:,:,:,8)));

    brainmask{n} = (im_diff{n} < 0.038) & abs(imallplus_all{n}(:,:,:,1)) > mean(mean(mean(abs(imallplus_all{n}(:,:,:,1)))))*.5;

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

figure
subplot(131)
imagesc(squeeze(abs(flip(im_masked{1}(:,:,46,1))))); colormap gray
subplot(132)
imagesc(imrotate(squeeze(abs(flip(im_masked{1}(:,50,:,1)))), 270));
subplot(133)
imagesc(imrotate(squeeze(abs(flip(im_masked{1}(50,:,:,1)))), 270));

%% setup data

B0 = 3;
phi_RF = .0562;
fit_thresh = [.011, .025];

% account for different z values if necessary or single slice test
for i = 1:3
    if exist('brainmask')
        im_masked{i} = im_masked{i}(:,:,:,1);
    else
        imallplus_all{i} = imallplus_all{i}(:,:,3:92,1); % remove top and bottom 2 slices for 3D UTE AFI cones
    end
end

n=2; % 2
    if exist('brainmask')
        im_uteplus_sample = im_masked{n}(:,:,:,1);            % brainmask passed into fitting
    else
        im_uteplus_sample = imallplus_all{n}(:,:,3:92,1);
    end
    
    fit_mask= repmat(abs(im_uteplus_sample) > fit_thresh(1), [1,1,1,size(imallplus_all{n},4)]);
for n=1:3 %3
    imallplus_all{n} = imallplus_all{n} .* fit_mask;
    im_uteplus{n} = imallplus_all{n}(:,:,:,1) .* fit_mask(:,:,:,1);
end



        imsize = size(im_uteplus_sample);
        
        Nechoes = size(imallplus_all{1}, 4);
        
        % only fit non-zero (can be used for masking) values
        I = find(im_uteplus{1});
        
        NI = length(I);
        
%       
        Snormall = max([im_uteplus{1}(:);im_uteplus{2}(:);im_uteplus{3}(:)]);
         % Snormall = max([im_uteplus{1}(:)]);
%         

imfit = cell(1,3); % cell(1, 3)
        imfitt = cell(1,3); % cell(1, 3)
for n=1:3 % 3
        imfit{n} = reshape(imallplus_all{n}, [prod(imsize) Nechoes]) / Snormall;
        
       imfitt{n} = imfit{n}(I,:);
end           

% for n = 1:3 % 3
%     TEin{n} = TE*1e-3; % s
% 
% end

% check for non-zero values and generate non-zero AFI volume

if exist('FA_map_ute_orientation')
    if size(find(FA_map_ute_orientation(I)), 1) < size(I, 1)
        x = FA_map_ute_orientation;%_ute_orientation;
        x(x <= 0) = AFI_nominal_flip;
        x(isnan(x)) = AFI_nominal_flip;
        nzFA_map = x;        
    end
    nzFA_map = FA_map_ute_orientation;
    imagesc(flip(nzFA_map(:,:,46))); colorbar % visualize non-zero FA map
end


%% fitting_multiscan_test

clear('fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3');
% make sure subet of voxels for fitting is all the same 
parfor Ix = 1: NI
% for Ix = [1 308701]
% disp( ['Ix: ', num2str(Ix) ] )
    
    Sin = cell(1,1); % cell(1, 3)
    for n = 1:3 %3
      Sin{n} = imfitt{n}(Ix,:); 
    end
    
    [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
        utebrain_multiscan_fitting_function(TEin, Sin, B0, phi_RF, 0);
       
    
    
%        disp([int2str(Ix) ' of ' int2str(NI)])
%          dbstop if error
end


savefname = '2020-02-03_multiscan-fitting_plus_dc';
writedata_flag = 1;

if writedata_flag
    save(savefname,'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
end   




%% fitting_T1_test

clear('fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3');


    

% make sure subet of voxels for fitting is all the same 
parfor Ix = 1: NI
% for Ix = [1 212965]
% disp( ['Ix: ', num2str(Ix) ] )
    
    Sin = cell(1,3);
    for n = 1:3
      Sin{n} = imfitt{n}(Ix,:); 
    end
    
    
    %B1_scale_map  = nzFA_map(I(Ix)) ./ AFI_nominal_flip;
    %aflips = flips .* B1_scale_map;
    
    [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
        utebrain_t1_fitting_function(TEin, Sin, flips, TR, B0, phi_RF, 1);
    
    
       %disp([int2str(Ix) ' of ' int2str(NI)])
         % dbstop if error
end


savefname = [scan_date '_T1-fitting_plus_aps_scale'];
writedata_flag = 1;

if writedata_flag
    save(savefname,'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
end

% generate fit maps
fit_filename = [datapath scan_date dateFolder filesep savefname];
fit_plus_filename = [datapath scan_date dateFolder filesep savefname]; 
AIC_thresh = 20;

fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh);

return
%% single  /ROI extraction
Itest = [24 39 43]; %


Iavg = 0;% [-1:1];

for n = 1:3
Sin{n} = squeeze(mean(mean(mean(imall_all{n}(Itest(1)+Iavg, Itest(2)+Iavg, Itest(3)+Iavg, :),1),2),3));
%S0plus = squeeze(mean(mean(mean(imallplus(Itest(1)+Iavg, Itest(2)+Iavg, Itest(3)+Iavg, :),1),2),3));

TEin{n} = TE*1e-3; % ms
% 
% % figure(98)
% % mpplot(TE, S0)
% % figure(99)
% % mpplot(TE, S0plus)
% 
end

%% This should work

 % keeps df, T2 fixed for a given component across experiments
[fit_result1, AIC1, fit_result2, AIC2, fit_result2m, AIC2m, fit_result3, AIC3] = ...
        utebrain_multiscan_fitting_function(TEin, Sin,3, 0, 1);
%     
%     %% combined T1 fitting - working now?
%     
     % keeps df, T2, T1 fixed for a given component across experiments
[fit_result1, AIC1, fit_result2, AIC2, fit_result2m, AIC2m, fit_result3, AIC3] = ...
        utebrain_t1_fitting_function(TEin, Sin,flips, TR, 3, 0, 1);
   
    
    
    %%
subplot(121)
imagesc(flip(fit_maps_0826_T1.T1(:,:,46))); colorbar

subplot(122)
imagesc(flip(fit_maps_0p4.T1(:,:,46))); colorbar

% subplot(122)
% histogram(fit_maps_0826_T1.T1)

% subplot(224)
% histogram(fit_maps_0p4.T1)
    