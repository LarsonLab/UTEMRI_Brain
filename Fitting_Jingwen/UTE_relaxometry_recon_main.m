clear; clc;

%% Add functions to path, set flags

addpath(genpath('/home/jyao3/020_UTE_Brain/UTEMRI_Brain'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/data/larson2/brain_uT2/orchestra-sdk-1.7-1.matlab'));
addpath(genpath('/home/plarson/matlab/3DUTE-recon'));

flag_plot = 0;
set(0,'DefaultAxesFOntSize',14);
set(0,'DefaultLineLineWidth',1.5);

%% specify raw data

subjectPath = '/data/larson2/brain_uT2/2020-06-02_3T_8TE';
cd(subjectPath);

pfile = dir(subjectPath);
pfile([pfile.isdir]) = [];
pfile(~cellfun(@isempty,strfind({pfile.name},'.mat'))) = [];
pfile(cellfun(@isempty,strfind({pfile.name},'P'))) = [];
pfile(~cellfun(@isempty,strfind({pfile.name},'recon'))) = [];

pfilename = {pfile(1).name, pfile(2).name, pfile(3).name};
TE_set = {'multi_utes-nd_8TE_18-deg.dat', 'multi_utes-nd_8TE_12-deg.dat', 'multi_utes-nd_8TE_06-deg.dat'};

%% reconstruct images

% 300 Hz offset
offset_frequency = 300;

% set up parameters
coils = 1:32;
undersamp = [];
skip = [];
freq_shift = offset_frequency;
echoes = [];
reg_coe = '-S -l2 -r0.01';
skip_calib_coil = 1;
cc_coil = 0;
rNecho = 0;
ind_echo_recon = 1;
espirit_recon = 0; % 1 for the older version

for ii = 2:length(pfilename)
    TE = get_TE(TE_set{ii});
    pfile = pfilename{ii};
    [im,header] = precon_3dute_pfile_bartv300_allec...
        (pfile, coils, undersamp, skip, freq_shift, echoes, ...
        reg_coe, skip_calib_coil, cc_coil, rNecho, ind_echo_recon, espirit_recon);
    
    % coil combination with COMPOSER
    phase = angle(im);
    magni = abs(im);
    phase_demod = phase - repmat(phase(:,:,:,:,1),[1 1 1 1 size(phase,5)]);
    compl_demod = magni.*exp(1i*phase_demod);
    compl_sum = squeeze(sum(compl_demod,4));
    
    save([pfilename{ii} '-NUFFT.mat'],'-v7.3',...
        'im','compl_sum','header','TE','offset_frequency');
end

% [im,header] = precon_3dute_pfile_bartv300_allec...
%     (pfile, coils, undersamp, skip, freq_shift, echoes, ...
%     reg_coe, skip_calib_coil, cc_coil, rNecho, ind_echo_recon, espirit_recon);
% save([pfilename{IndData} '-SENSE.mat'],'-v7.3','im','header','TE','offset_frequency');

%% plot COMPOSER combined phase

if flag_plot
    
    im_slice1 = compl_sum(:,:,30,:);
    im_slice2 = compl_sum(:,:,46,:);
    im_slice3 = compl_sum(:,:,60,:);
    
    figure('Position',[100,100,1300,900]);
    subplot(311)
    imagesc(angle(im_slice1(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 30');
    subplot(312)
    imagesc(angle(im_slice2(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 46');
    subplot(313)
    imagesc(angle(im_slice3(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 60');
    export_fig(sprintf('COMPOSER_phase'),'-png','-transparent'); close;
    
end

%% espirit combination

if flag_plot
    
    load([pfilename{ii} '-csreconallec_l2_r0p01.mat']);
    
    im_slice1 = imallplus(:,:,30,:);
    im_slice2 = imallplus(:,:,46,:);
    im_slice3 = imallplus(:,:,60,:);
    
    figure('Position',[100,100,1300,900]);
    subplot(311)
    imagesc(angle(im_slice1(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 30');
    subplot(312)
    imagesc(angle(im_slice2(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 46');
    subplot(313)
    imagesc(angle(im_slice3(:,:)),[-pi pi]); axis equal; axis tight; axis off
    title('Slice 60');
    export_fig(sprintf('Coil_comb_phase'),'-png','-transparent'); close;
    
end
