clear; clc;

%% Add functions to path, set flags

addpath(genpath('/home/jyao3/020_UTE_Brain/UTEMRI_Brain'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/data/larson2/brain_uT2/orchestra-sdk-1.7-1.matlab'));

flag_writefiles = 0;
flag_testOneSlice = 1; slice_number = 46;

plot_flag = 0;
set(0,'DefaultAxesFOntSize',14);
set(0,'DefaultLineLineWidth',1.5);

%% Set parameters

B0 = 3; % Tesla
phi_RF = .0562;
TR = 9.1e-3; % s
methylene_freq_est = 3.5 * B0 * 42.57e-3; % kHz

%% Set input file path

datapath = '/data/larson2/brain_uT2/';
scan_date = '2020-06-02';
suffix = '-NUFFT.mat';
dateFolder = '_3T_8TE/';

Pfile_18deg = [datapath scan_date dateFolder '/P69632.7_06021135'];
Pfile_12deg = [datapath scan_date dateFolder '/P70656.7_06021148'];
Pfile_06deg = [datapath scan_date dateFolder '/P73216.7_06021202'];
Pfile_AFI = [datapath scan_date dateFolder '/P74240.7_06021210'];

%% Autoprescan correction data
% Needed if autoprescan is done instead of manual prescan

% read in auto prescan values to correct scaling
pfile_18 = GERecon('Pfile.Load', Pfile_18deg);
header_18 = GERecon('Pfile.Header', pfile_18);

pfile_12 = GERecon('Pfile.Load', Pfile_12deg);
header_12 = GERecon('Pfile.Header', pfile_12);

pfile_06 = GERecon('Pfile.Load', Pfile_06deg);
header_06 = GERecon('Pfile.Header', pfile_06);

auto_prescan_values = struct('aps_r1', {}, 'aps_r2', {}, 'aps_tg', {});
auto_prescan_values(1).aps_r1 = ...
    [header_18.PrescanHeader.aps_r1, header_12.PrescanHeader.aps_r1, header_06.PrescanHeader.aps_r1];
auto_prescan_values(1).aps_r2 = ...
    [header_18.PrescanHeader.aps_r2, header_12.PrescanHeader.aps_r2, header_06.PrescanHeader.aps_r2];
auto_prescan_values(1).aps_tg = ...
    [header_18.PrescanHeader.aps_tg, header_12.PrescanHeader.aps_tg, header_06.PrescanHeader.aps_tg];

disp(['  -- Aotoprescan values r1: ' num2str(auto_prescan_values.aps_r1)]);
disp(['  -- Aotoprescan values r2: ' num2str(auto_prescan_values.aps_r2)]);
disp(['  -- Aotoprescan values tg: ' num2str(auto_prescan_values.aps_tg)]);

%% Load data

% FA1: 18-deg
load([Pfile_18deg suffix],'compl_sum','TE');

imallplus_all{3} = compl_sum;
TEin{3} = TE;
flips(3) = 18*pi/180;

% FA2: 12-deg
load([Pfile_12deg suffix],'compl_sum','TE');

imallplus_all{2} = compl_sum .* (10 ^ (3 * (header_18.PrescanHeader.aps_r1 - header_12.PrescanHeader.aps_r1) / 20)) ...
    .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_12.PrescanHeader.aps_r2));% * 1e4;
flips(2) = 12*pi/180;
TEin{2} = TE;

% FA3: 6-deg
load([Pfile_06deg suffix],'compl_sum','TE');

imallplus_all{1} = compl_sum .* 10 ^(3 * (header_18.PrescanHeader.aps_r1 - header_06.PrescanHeader.aps_r1) /20) ...
    .* (2 ^ (header_18.PrescanHeader.aps_r2 - header_06.PrescanHeader.aps_r2));% * 1e4;
flips(1) = 6*pi/180;
TEin{1} = TE;

for j = 1:3
    TEin{j} = TEin{j} * 1e-3;
end

disp(['  -- FA  6 TE (ms): ' num2str(TEin{1})]);
disp(['  -- FA 12 TE (ms): ' num2str(TEin{2})]);
disp(['  -- FA 18 TE (ms): ' num2str(TEin{3})]);

%% Visualize flip angle volumes and check scaling

disp(['Size of each FA volume: ' num2str(size(imallplus_all{1}))]);

figure('Position',[100 100 1200 300]);
subplot(131)
imagesc(flip(abs(imallplus_all{3}(:,:,46,1)))); colormap gray; colorbar
subplot(132)
imagesc(flip(abs(imallplus_all{2}(:,:,46,1))));colormap gray; colorbar
subplot(133)
imagesc(flip(abs(imallplus_all{1}(:,:,46,1))));colormap gray; colorbar

%% B1 correction

% AFI_nominal_flip = 44 * pi/180; % UTE AFI radial
AFI_nominal_flip = 45 * pi/180; % UTE AFI cones

% load in UTE AFI volumes
load([Pfile_AFI '_UTE_AFI.mat']);
AFI_map_1 = squeeze(recon_grid(:,:,1,:));
AFI_map_2 = squeeze(recon_grid(:,:,2,:));

% generate flip angle map - AFI
S1 = abs(AFI_map_1); S2 = abs(AFI_map_2); % magnitude images
TR1 = 7; TR2 = 35; % ms
n = TR2 / TR1; r = S2 ./ S1;
FA_map = abs(acos((n * r - 1) ./ (n - r))); % rad

disp(['Size of AFI volume: ' num2str(size(FA_map))]);

FA_map_resize = FA_map(:, 6:91, :); % make sure dimensions match UTE image dimensions
% FA_map_deg = FA_map_resize * 180/pi; % convert to degrees if necessary
% FA_map_scaled = FA_map_resize ./ AFI_nominal_flip;

% Gaussian smoothing
FA_map_resize = imgaussfilt3(FA_map_resize,1);

figure('Position',[100 100 800 300]);
subplot(121)
imagesc(FA_map_resize(:,:,46),[0.5 1.2]); colormap gray; colorbar
subplot(122)
histogram(FA_map_resize,[0:0.02:2]); hold on;
scatter(AFI_nominal_flip,0,'*');

% hard code correction test
FA_map_resize_cor = FA_map_resize + 0.35;
% figure(2)
% histogram(FA_map_resize_cor)

%% AFI B1 map orientation correction to match UTE relaxometry seq
% FA_map_ute_orientation = flip(flip(imrotate(FA_map_resize_cor, 270), 2), 3); % 3D UTE AFI radial
% FA_map_ute_orientation = flip(FA_map_resize, 3);
FA_map_ute_orientation = flip(flip(FA_map_resize, 3), 1); % 3D UTE AFI cones

% check orientation
figure('Position',[100 100 800 500])
subplot(231)
imagesc(abs(flip(imallplus_all{1}(:,:,46)))); colormap gray;
subplot(232)
imagesc(imrotate(abs(squeeze(imallplus_all{1}(:,50,:,1))), 270)); colormap gray
subplot(233)
imagesc(imrotate(abs(squeeze(imallplus_all{1}(50,:,:,1))), 270)); colormap gray

subplot(234)
imagesc(flip(FA_map_ute_orientation(:,:,46)), [0 1.5]); colormap gray;
subplot(235)
imagesc(imrotate(squeeze(FA_map_ute_orientation(:,50,:,1)), 270), [0 1.5]); colormap gray
subplot(236)
imagesc(imrotate(squeeze(FA_map_ute_orientation(50,:,:,1)), 270), [0 1.5]); colormap gray

%% brainmask

% Jingwen - add bet2
addpath('/netopt/share/ese/ESE_DV26.0_R01/tools/matlab/read_MR/');

brainmask = cell(1, 3);

% read in dimension info
FOV = readFOVfromPfile(Pfile_18deg);
MatSize = size(imallplus_all{1}(:,:,:,1));
voxSize = FOV./MatSize; % not 2x2x2?

for n = 1:3
    
    brainmask{n} = BET(abs(imallplus_all{n}(:,:,:,end)),MatSize,voxSize);
    
    %     intensity based thresholding
    %     brainmask{n} = abs(imallplus_all{n}(:,:,:,1)) ...
    %         > mean(mean(mean(abs(imallplus_all{n}(:,:,:,1)))))*.5;
    %         % & (im_diff{n} < 0.038);
    %
    %     se = strel('disk', 5);
    %     brainmask{n} = imopen(brainmask{n}, se);
    %
    %     brainmask{n}(54:end,:,60:end) = 0;
    %     brainmask{n}(73:end,:,55:end) = 0;
    %     brainmask{n}(95:end,:,44:end) = 0;
    %     brainmask{n}(:,:,73:end) = 0;
    %     brainmask{n}(:,:,1:10) = 0;
    
end

brainmask_comb = brainmask{1} & brainmask{2} & brainmask{3};

im_masked{1} = imallplus_all{3}(:,:,:,1) .* brainmask_comb;

figure('Position',[100 100 1200 600])
subplot(231)
imagesc(squeeze(abs(flip(brainmask_comb(:,:,46,1))))); colormap gray
subplot(232)
imagesc(imrotate(squeeze(abs(flip(brainmask_comb(:,50,:,1)))), 270));
subplot(233)
imagesc(imrotate(squeeze(abs(flip(brainmask_comb(50,:,:,1)))), 270));
subplot(234)
imagesc(squeeze(abs(flip(im_masked{1}(:,:,46,1))))); colormap gray
subplot(235)
imagesc(imrotate(squeeze(abs(flip(im_masked{1}(:,50,:,1)))), 270));
subplot(236)
imagesc(imrotate(squeeze(abs(flip(im_masked{1}(50,:,:,1)))), 270));

export_fig('brain_mask','-png','-transparent'); close;

%% phase unwrapping

im_masked = cell(1, 3);
im_masked_magni = cell(1, 3);
im_masked_phase = cell(1, 3);
Nechoes = size(imallplus_all{1},4);

% mask data
brainmask_cut = brainmask_comb(:,:,3:92); % match slices of AFI data
for ii = 1:3
    im_masked{ii} = imallplus_all{ii}(:,:,3:92,:);
    im_masked{ii} = im_masked{ii}.* repmat(brainmask_cut,[1 1 1 Nechoes]);
    im_masked_magni{ii} = abs(im_masked{ii});
    Inputs.Mask = double(brainmask_cut);
    for ee = 1:Nechoes
        Inputs.Phase = angle(im_masked{ii}(:,:,:,ee));
        Unwrapped = SEGUE(Inputs);
        im_masked_phase{ii}(:,:,:,ee) = Unwrapped;
    end
    im_masked{ii} = im_masked_magni{ii}.*exp(1i*im_masked_phase{ii});
end

figure('Position',[100 100 1200 600])
subplot(231)
imagesc(squeeze(angle(flip(im_masked{1}(:,:,46,end))))); colormap gray
subplot(232)
imagesc(imrotate(squeeze(angle(flip(im_masked{1}(:,50,:,end)))), 270));
subplot(233)
imagesc(imrotate(squeeze(angle(flip(im_masked{1}(50,:,:,end)))), 270));
subplot(234)
imagesc(squeeze((flip(im_masked_phase{1}(:,:,46,end))))); colormap gray
subplot(235)
imagesc(imrotate(squeeze((flip(im_masked_phase{1}(:,50,:,end)))), 270));
subplot(236)
imagesc(imrotate(squeeze((flip(im_masked_phase{1}(50,:,:,end)))), 270));

export_fig('image_phase','-png','-transparent'); close;

%% Setup data for fitting

% fit_thresh = [.011 .025];

% normalize data
Smax = max(abs([im_masked{1}(:);im_masked{2}(:);im_masked{3}(:)]));
im_norm = cell(1, 3);
for ii = 1:3
    im_norm{ii} = im_masked{ii}/Smax;
end

% select slice or slab
if ~flag_testOneSlice
    slice_number = 1:size(brainmask_cut,3);
end
brainmask_cut_test = brainmask_cut(:,:,slice_number);
im_norm_test = cell(1,3);
for i = 1:3
    im_norm_test{i} = im_norm{i}(:,:,slice_number,:);
end
FA_map_ute_test = FA_map_ute_orientation(:,:,slice_number);

% only fit within the mask
I = find(brainmask_cut_test);
imfit = cell(1, 3);
for i = 1:3
    imfit{i} = reshape(im_norm_test{i},[],Nechoes);
    imfit{i} = imfit{i}(I,:);
end

disp([' -- Fitting voxel number: ' num2str(length(I))]);

% check FAmap values and fill in missing ones
I_FA = find(brainmask_cut_test & (isnan(FA_map_ute_test) | FA_map_ute_test == 0));
if ~isempty(I_FA)
    FA_map_ute_test(I_FA) = AFI_nominal_flip;
end
FAmapfit = FA_map_ute_test(I);

%% Select ROI

% figure; imagesc((abs(im_norm_test{3}(:,:,end))));
% h = imellipse;
% ROI_WM = h.createMask;
% close;
% 
% figure; imagesc((abs(im_norm_test{3}(:,:,end))));
% h = imellipse;
% ROI_GM = h.createMask;
% close;
% 
% save('test_roi.mat','ROI_WM','ROI_GM');

%% fitting 1 component

% clear('fitting_result');
% 
% % loop through voxels for fitting
% for Ix = 1000 % 1:length(I)
%     
%     if mod(Ix,100) == 0
%         disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
%     end
%     
%     Sin = cell(1,3);
%     for i = 1:3
%         Sin{i} = imfit{i}(Ix,:);
%     end
%     
%     if plot_flag
%         figure('Position',[100 100 400 500]);
%         plot(TEin{1},abs(Sin{1})); hold on;
%         plot(TEin{2},abs(Sin{2}));
%         plot(TEin{3},abs(Sin{3}));
%         legend({'FA = 6deg','12deg','18deg'},'box','off','location','best');
%         export_fig('fitplot_signal','-png','-transparent'); close;
%     end
%     
%     B1_scale  = FAmapfit(Ix) ./ AFI_nominal_flip;
%     aflips = flips .* B1_scale;
%     
%     [fitting_result_pre(Ix)] = ...
%         UTE_fitting_function_step1(TEin, Sin, aflips, TR, B0, phi_RF, plot_flag);
%     
% end
% 
% %% phase unwrapping
% 
% fit_maps = UTE_save_maps(fitting_result_pre, 'comp1T1',brainmask_cut_test, I);
% phase1 = fit_maps{1}.phi1; phase{1}.unwrap = unwrap_phase(phase1);
% phase2 = fit_maps{1}.phi2; phase{2}.unwrap = unwrap_phase(phase2);
% phase3 = fit_maps{1}.phi3; phase{3}.unwrap = unwrap_phase(phase3);
% 
% figure('position',[100 100 1600 600]); 
% clear ax;
% ax(1) = subplot(241); 
% imagesc(flip(phase1),[-pi pi]); colorbar; title('Comp1 phi1');
% ax(2) = subplot(242); 
% imagesc(flip(phase2),[-pi pi]); colorbar; title('Comp1 phi2');
% ax(3) = subplot(243); 
% imagesc(flip(phase3),[-pi pi]); colorbar; title('Comp1 phi3');
% ax(4) = subplot(244); 
% phiSD = nanstd(cat(3,phase1,phase2,phase3),[],3);
% imagesc(flip(phiSD),[-pi pi]/10); colorbar; title('Comp1 phi SD (x10)');
% ax(5) = subplot(245); 
% imagesc(flip(phase{1}.unwrap),[-pi pi]); colorbar; title('Comp1 phi1');
% ax(6) = subplot(246); 
% imagesc(flip(phase{2}.unwrap),[-pi pi]); colorbar; title('Comp1 phi2');
% ax(7) = subplot(247); 
% imagesc(flip(phase{3}.unwrap),[-pi pi]); colorbar; title('Comp1 phi3');
% ax(8) = subplot(248); 
% phiSDuw = nanstd(cat(3,phase{1}.unwrap,phase{2}.unwrap,phase{3}.unwrap),[],3);
% imagesc(flip(phiSDuw),[-pi pi]/10); colorbar; title('Comp1 phi SD (x10)');
% axis(ax,'off','equal','tight');
% 
% export_fig('phase_unwrap','-png','-transparent'); close;
% 
% phi_est = cell(1, 3);
% for i = 1:3
%     phi_est{i} = phase{i}.unwrap(I);
% end

%% fitting multiscan T2* and T1

% save('test_IDEAL_slice.mat','I','imfit','TEin','FAmapfit',...
%     'AFI_nominal_flip','flips','brainmask_cut_test','im_norm_test');

clear('fitting_result');

plot_flag = 1;

general_opts.TR = TR;
general_opts.B0 = B0;
general_opts.phi_RF = phi_RF;
general_opts.fixPhi = 0;
general_opts.singlePhi = 1;
general_opts.fixDf = 0;
general_opts.fixT1comp2 = 0; general_opts.T1median = 0.5259;
general_opts.plot_flag = plot_flag;
general_opts.methylene_freq_est = 3.5*B0*42.57e-3; % kHz % 3.5 2.31

% loop through voxels for fitting
for Ix = 500 % 1:length(I)
    
    if mod(Ix,100) == 0
        disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
    end
    
    Sin = cell(1,3);
    for i = 1:3
        Sin{i} = imfit{i}(Ix,:);
    end
    
    if plot_flag
        figure('Position',[100 100 400 500]);
        plot(TEin{1},abs(Sin{1})); hold on;
        plot(TEin{2},abs(Sin{2}));
        plot(TEin{3},abs(Sin{3}));
        legend({'FA = 6deg','12deg','18deg'},'box','off','location','best');
        export_fig('fitplot_signal','-png','-transparent'); close;
    end
    
    B1_scale  = FAmapfit(Ix) ./ AFI_nominal_flip;
    aflips = flips .* B1_scale; 
    
    % without phase unwrapping
    [fitting_result(Ix)] = ...
        UTE_joint_fitting_function(TEin, Sin, aflips, general_opts);
    
    % from Github script
%     [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), ~, ~, ~, ~] = ...
%         utebrain_t1_fitting_function(TEin, Sin, aflips, TR, B0, phi_RF, 0);
    
    % with phase unwrapping
%     PHIin = zeros(1,3);
%     for i = 1:3
%         PHIin(i) = phi_est{i}(Ix);
%     end
%     
%     [fitting_result(Ix)] = ...
%         UTE_fitting_function_step2(TEin, Sin, aflips, general_opts, fitting_result_pre(Ix), PHIin);

end

% save('test_IDEAL_data.mat','TEin','Sin','aflips');

%% save and plot parameter maps

AIC_thresh = 20;
fit_maps = UTE_save_maps(fitting_result, 'comp2T1',...
    brainmask_cut_test, I, AIC_thresh);

% parameter maps
figure('position',[100 100 1600 600]); 
ax(1) = subplot(241);
imagesc(flip(fit_maps{1}.rho),[0 20]); colorbar; title('Comp1 rho');
ax(2) = subplot(242); 
imagesc(flip(fit_maps{2}.rho),[0 2]); colorbar; title('Comp2 rho');
ax(3) = subplot(243); 
imagesc(flip(fit_maps{2}.rho_frac),[0 0.3]); colorbar; title('Comp2 fraction');
ax(4) = subplot(244); 
imagesc(flip(fit_maps{1}.AIC)); colorbar; title('AIC');
ax(5) = subplot(245); 
imagesc(flip(fit_maps{1}.T1),[0 5]); colorbar; title('Comp1 T1');
ax(6) = subplot(246); 
imagesc(flip(fit_maps{2}.T1),[0 1]); colorbar; title('Comp2 T1');
ax(7) = subplot(247); 
imagesc(flip(fit_maps{1}.T2),[0 100]); colorbar; title('Comp1 T2');
ax(8) = subplot(248); 
imagesc(flip(fit_maps{2}.T2),[0 1]); colorbar; title('Comp2 T2');
axis(ax,'off','equal','tight');

% CV for parameter estimations
figure('position',[100 100 1600 600]); 
ax(1) = subplot(241); 
imagesc(flip(fit_maps{1}.rho_cv),[0 20]); colorbar; title('Comp1 rho CV');
ax(2) = subplot(242); 
imagesc(flip(fit_maps{2}.rho_cv),[0 50]); colorbar; title('Comp2 rho CV');
ax(3) = subplot(243); 
imagesc(flip(fit_maps{1}.T1_cv),[0 20]); colorbar; title('Comp1 T1 CV');
ax(4) = subplot(244); 
imagesc(flip(fit_maps{2}.T1_cv),[0 70]); colorbar; title('Comp2 T1 CV');
ax(5) = subplot(245); 
imagesc(flip(fit_maps{1}.T2_cv),[0 30]); colorbar; title('Comp1 T2 CV');
ax(6) = subplot(246); 
imagesc(flip(fit_maps{2}.T2_cv),[0 50]); colorbar; title('Comp2 T2 CV');
ax(7) = subplot(247); 
imagesc(flip(fit_maps{2}.df/B0/42.57e-3),[0 5]); colorbar; title('Comp2 df');
ax(8) = subplot(248); 
imagesc(flip(fit_maps{2}.df_err/B0/42.57e-3),[0 5]); colorbar; title('Comp2 df err');
axis(ax,'off','equal','tight');

% phi maps
clear ax
figure('position',[100 100 1600 600]); 
ax(1) = subplot(241); 
imagesc(flip(fit_maps{1}.phi1),[-pi pi]); colorbar; title('Comp1 phi1');
ax(2) = subplot(242); 
imagesc(flip(fit_maps{1}.phi2),[-pi pi]); colorbar; title('Comp1 phi2');
ax(3) = subplot(243); 
imagesc(flip(fit_maps{1}.phi3),[-pi pi]); colorbar; title('Comp1 phi3');
ax(4) = subplot(244); 
phiSD = nanstd(cat(3,fit_maps{1}.phi1, fit_maps{1}.phi2, fit_maps{1}.phi3),[],3);
imagesc(flip(phiSD),[-pi pi]/10); colorbar; title('Comp1 phi SD (x10)');
ax(5) = subplot(245); 
imagesc(flip(fit_maps{2}.phi1),[-pi pi]); colorbar; title('Comp2 phi1');
ax(6) = subplot(246); 
imagesc(flip(fit_maps{2}.phi2),[-pi pi]); colorbar; title('Comp2 phi2');
ax(7) = subplot(247); 
imagesc(flip(fit_maps{2}.phi3),[-pi pi]); colorbar; title('Comp2 phi3');
ax(8) = subplot(248); 
phiSD = nanstd(cat(3,fit_maps{2}.phi1, fit_maps{2}.phi2, fit_maps{2}.phi3),[],3);
imagesc(flip(phiSD),[-pi pi]/10); colorbar; title('Comp2 phi SD (x10)');
axis(ax,'off','equal','tight');

export_fig('fitmap3','-png','-transparent'); close;
export_fig('fitmap2','-png','-transparent'); close;
export_fig('fitmap1','-png','-transparent'); close;

% save('test_new.mat','fitting_result','fit_maps');

%% plot df stats

figure('position',[100 100 300 300]); 
histogram(nonzeros(fit_maps{2}.df/B0/42.57e-3),-2:0.3:6); hold on;
xlabel('df (ppm)'); ylabel('Voxel count');
medDf = nanmedian(nonzeros(fit_maps{2}.df/B0/42.57e-3));
scatter(medDf,0,'s','filled');
title(sprintf('Median df %.2f ppm',medDf));
export_fig('df_hist','-png','-transparent'); close;

%% plot T1 stats

figure('position',[100 100 300 300]); 
scatter(fit_maps{2}.T1(I),fit_maps{2}.T1_cv(I),'.'); hold on;
xlabel('Comp2 T1'); ylabel('Comp1 T1 CV (%)');
ylim([0 100]); xlim([0 1]);

thresh = prctile(fit_maps{2}.T1_cv(I),10);
line([0 1], [thresh thresh], 'Color', 'k');

T1belowTh = fit_maps{2}.T1(fit_maps{2}.T1_cv < thresh);
medianT1 = nanmedian(T1belowTh);
scatter(medianT1, thresh, 'rd', 'filled');

%% plot ROI stats

% load('test_roi.mat');

% figure('position',[100 100 600 600]); 
% h = imagesc(flip((abs(im_norm_test{3}(:,:,end))))); hold on;
% h2 = imagesc(flip(ones(size(im_norm_test{3}(:,:,end)))));
% axis tight; axis equal; axis off;
% set(h2, 'AlphaData', flip(ROI_WM+ROI_GM)*0.5);

comp = 1; fd = 'rho_cv';
map = fit_maps{comp}.(fd);

WMval = nonzeros(map.*ROI_WM); WMval(isnan(WMval)) = [];
GMval = nonzeros(map.*ROI_GM); GMval(isnan(GMval)) = [];
histBin = linspace(min([WMval; GMval]), max([[WMval; GMval]]),10);

figure('position',[100 100 300 300]); 
histogram(WMval,histBin); hold on;
histogram(GMval,histBin);
legend({'WM','GM'},'box','off');
xlabel(sprintf('Component %i %s', comp, fd), 'Interpreter', 'none');
export_fig(sprintf('comp%i%s',comp,fd),'-png','-transparent'); close;
