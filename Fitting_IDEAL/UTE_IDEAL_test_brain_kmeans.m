clc; clear;
warning('off');

%% add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/home/jyao3/020_UTE_Brain'));

%% load data

load('test_IDEAL_slice.mat','TEin','flips');
dataPath = '/working/larson/jingwen/030_UTE_Brain/test_data/';

nii = load_nii([dataPath '/Img_magni.nii.gz']);
Mag = double(nii.img);
nii = load_nii([dataPath '/Img_phase.nii.gz']);
Pha = double(nii.img);
nii = load_nii([dataPath '/brain_mask.nii.gz']);
mask = double(nii.img);
nii = load_nii([dataPath '/B1map.nii.gz']);
B1map = double(nii.img).*mask;

Img = Mag.*exp(1i*Pha).*repmat(mask,[1 1 1 size(Mag,4)]);

%% select slice and plot

sliceN = 46;

Img_slice = squeeze(Img(:,:,sliceN,:));
mask_slice = mask(:,:,sliceN);
B1map_slice = B1map(:,:,sliceN);

figure('Position',[100 100 400 400]);
subplot(221)
imagesc(flip(abs(Img_slice(:,:,3)))); axis off; axis tight; axis equal;
subplot(222)
imagesc(flip(angle(Img_slice(:,:,3)))); axis off; axis tight; axis equal;
subplot(223)
imagesc(flip(abs(B1map_slice))); axis off; axis tight; axis equal;
subplot(224)
imagesc(flip(abs(mask_slice))); axis off; axis tight; axis equal;
export_fig(sprintf('SliceImg'),'-png','-transparent'); close;

%% set up parameters

B0 = 3; % Tesla
TR = 9.1e-3; % s
methylene_freq_est = 3.5 * B0 * 42.57e-3; % kHz

% from 10.1073/PNAS.1115107109
multipeak.alpha = [0.743 0.124 0.111 0.021];
multipeak.alpha = multipeak.alpha/sum(multipeak.alpha);
multipeak.chemshift = [1.55 0.9 3.2 1.3]-4.8; % ppm
multipeak.chemshift = multipeak.chemshift * B0 * 42.57e-3;

%% prepare data

TEall = transpose([TEin{:}])/1000; % s
Necho = length(TEall);
Ncomp = 2;
Nacq = length(TEin);

Img_mask = reshape(Img_slice,[],Necho);
I = find(mask_slice > 0);
Sinput = Img_mask(I,:); % Nvoxel x 24

% normalization
Sinput = Sinput/max(abs(Sinput(:)));

B1input  = B1map_slice(I);
FAall = repmat(flips,[Necho/length(TEin),1]); FAall = FAall(:);
FAinput = repmat(FAall',[length(I),1]) .* repmat(B1input,[1,Necho]); % Nvoxel x 24

%% IDEAL

brainmask_cut_test = mask_slice>0;

EstValue.Rho = zeros(length(I),2);
EstValue.Rho_real = zeros(length(I),2);
EstValue.Rho_imag = zeros(length(I),2);
EstValue.R2s = zeros(length(I),2);
EstValue.R1 = zeros(length(I),2);
EstValue.Phi = zeros(length(I),1);
EstValue.Rho_frac = zeros(length(I),1);
EstValue.Rmse = zeros(length(I),1);
EstValue.Iter = zeros(length(I),1);

idx = kmeans([real(Sinput) imag(Sinput)],3);
% figure; imagesc(flip(img_from_fit(idx, brainmask_cut_test)));

for m = 1:5
    
    indN = 3^m;
    idx = kmeans([real(Sinput) imag(Sinput)],indN);
    fprintf('K-means iteration %i Nclusters %i \n',m,indN);
    
    Rho = zeros(indN,2);
    Rho_real = zeros(indN,2);
    Rho_imag = zeros(indN,2);
    R2s = zeros(indN,2);
    R1 = zeros(indN,2);
    Phi = zeros(indN,1);
    Rho_frac = zeros(indN,1);
    Rmse = zeros(indN,1);
    Iter = zeros(indN,1);

    initCond = cell(1,indN);
    parfor ii = 1:indN
        Smean = transpose(mean(Sinput(idx == ii,:),1));
        FAmean = transpose(mean(FAinput(idx == ii,:),1));
        
        if m == 1
            initCond{ii}.phi = 60; % Hz
            initCond{ii}.R2s = [30; 3000]; % s-1
            initCond{ii}.R1 = [0.5; 3]; % s-1
        else
            initCond{ii}.phi = nanmean(EstValue.Phi(idx == ii),1);
            initCond{ii}.R2s = nanmean(EstValue.R2s(idx == ii,:),1)';
            initCond{ii}.R1 = nanmean(EstValue.R1(idx == ii,:),1)';
            if sum(isnan(initCond{ii}.phi)) ~= 0
                initCond{ii}.phi = 60;
            end
            if sum(isnan(initCond{ii}.R2s)) ~= 0
                initCond{ii}.R2s = [30; 3000];
            end
            if sum(isnan(initCond{ii}.R1)) ~= 0
                initCond{ii}.R1 = [0.5; 3];
            end
        end
        
        % IDEAL iterations
        iterMax = 100;
        thres = 0.1;
        [fitparam, ~, rmse_Iter, S_fit, iter] = ...
            IDEALfunc(Smean, TEall, FAmean, TR, multipeak, initCond{ii}, iterMax, thres);
        
        % assign component according to R2
        [fitparam.R2s,indComp] = sort(fitparam.R2s);
        fitparam.R1 = fitparam.R1(indComp);
        
        % save estimation values into arrays
        rho_mag = fitparam.rho(1:2:end) + 1i*fitparam.rho(2:2:end);
        rho_mag = rho_mag(indComp);
        
        Rho(ii,:) = abs(rho_mag)';
        Rho_real(ii,:) = real(rho_mag)';
        Rho_imag(ii,:) = imag(rho_mag)';
        R2s(ii,:) = fitparam.R2s';
        R1(ii,:) = fitparam.R1';
        Phi(ii) = fitparam.phi;
        Rho_frac(ii) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));

        Rmse(ii) = rmse_Iter(end);
        Iter(ii) = iter;
    end
    
    for ii = 1:indN
        EstValue.Rho(idx == ii,:) = repmat(Rho(ii,:),[sum(idx == ii),1]);
        EstValue.Rho_real(idx == ii,:) = repmat(Rho_real(ii,:),[sum(idx == ii),1]);
        EstValue.Rho_imag(idx == ii,:) = repmat(Rho_imag(ii,:),[sum(idx == ii),1]);
        EstValue.R2s(idx == ii,:) = repmat(R2s(ii,:),[sum(idx == ii),1]);
        EstValue.R1(idx == ii,:) = repmat(R1(ii,:),[sum(idx == ii),1]);
        EstValue.Phi(idx == ii) = Phi(ii);
        EstValue.Rho_frac(idx == ii) = Rho_frac(ii);
        EstValue.Rmse(idx == ii) = Rmse(ii);
        EstValue.Iter(idx == ii) = Iter(ii);
    end
    
    plot_EstMap(EstValue, brainmask_cut_test, idx);
    export_fig(sprintf(['IDEAL_slice2_iter' num2str(m)]),'-png','-transparent'); close; pause(1);
    export_fig(sprintf(['IDEAL_slice_iter' num2str(m)]),'-png','-transparent'); close;
    
    % smooth filter the B0 map
    EstMap.Phi_map = img_from_fit(EstValue.Phi(:,1), brainmask_cut_test);
    EstMap.Phi_map(isnan(EstMap.Phi_map)) = 0;
    EstMap.Phi_map = imgaussfilt(EstMap.Phi_map,2);
    EstValue.Phi = EstMap.Phi_map(brainmask_cut_test);
    
    % avoid extreme values in T1 comp2 and T2 comp1
%     prc = prctile(EstValue.R1(:,2),[2.5 97.5]);
%     EstValue.R1(EstValue.R1(:,2) < prc(1),2) = prc(1);
%     EstValue.R1(EstValue.R1(:,2) > prc(2),2) = prc(2);
%     prc = prctile(EstValue.R2s(:,1),[2.5 97.5]);
%     EstValue.R2s(EstValue.R2s(:,1) < prc(1),2) = prc(1);
%     EstValue.R2s(EstValue.R2s(:,1) > prc(2),2) = prc(2);
    
end

%% Full resolution

Rho = zeros(length(I),2);
Rho_real = zeros(length(I),2);
Rho_imag = zeros(length(I),2);
R2s = zeros(length(I),2);
R1 = zeros(length(I),2);
Phi = zeros(length(I),1);
Rho_frac = zeros(length(I),1);
Rmse = zeros(length(I),1);
Iter = zeros(length(I),1);

initCond = cell(length(I),1);
parfor Ix = 1:length(I)
    
    if mod(Ix,100) == 0
        disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
    end
    
    S0all = transpose(Sinput(Ix,:));
    FAall = transpose(FAinput(Ix,:));
    
    initCond{Ix}.phi = EstValue.Phi(Ix);
    initCond{Ix}.R2s = EstValue.R2s(Ix,:)';
    initCond{Ix}.R1 = EstValue.R1(Ix,:)';
    if sum(isnan(initCond{Ix}.phi)) ~= 0
        initCond{Ix}.phi = 0;
    end
    if sum(isnan(initCond{Ix}.R2s)) ~= 0
        initCond{Ix}.R2s = [30; 3000];
    end
    if sum(isnan(initCond{Ix}.R1)) ~= 0
        initCond{Ix}.R1 = [0.5; 3];
    end
    
    % IDEAL iterations
    iterMax = 200;
    thres = 0.1;
    [fitparam, ~, rmse_Iter, S_fit, iter] = ...
        IDEALfunc(S0all, TEall, FAall, TR, multipeak, initCond{Ix}, iterMax, thres);
    
    % assign component according to R2
    [fitparam.R2s,indComp] = sort(fitparam.R2s);
    fitparam.R1 = fitparam.R1(indComp);
    
    % save estimation values into arrays
    rho_mag = fitparam.rho(1:2:end) + 1i*fitparam.rho(2:2:end);
    rho_mag = rho_mag(indComp);
    
    Rho(Ix,:) = abs(rho_mag)';
    Rho_real(Ix,:) = real(rho_mag)';
    Rho_imag(Ix,:) = imag(rho_mag)';
    R2s(Ix,:) = fitparam.R2s';
    R1(Ix,:) = fitparam.R1';
    Phi(Ix) = fitparam.phi;
    Rho_frac(Ix) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));
    
    Rmse(Ix) = rmse_Iter(end);
    Iter(Ix) = iter;
    
end

FinalValue.Rho = Rho;
FinalValue.Rho_real = Rho_real;
FinalValue.Rho_imag = Rho_imag;
FinalValue.R2s = R2s;
FinalValue.R1 = R1;
FinalValue.Phi = Phi;
FinalValue.Rho_frac = Rho_frac;
FinalValue.Rmse = Rmse;
FinalValue.Iter = Iter;

plot_EstMap(FinalValue, brainmask_cut_test, []);
export_fig(sprintf('IDEAL_slice2_full'),'-png','-transparent'); close; pause(1);
export_fig(sprintf('IDEAL_slice_full'),'-png','-transparent'); close;

%%

EstMapReal.Rho1_map = flip(img_from_fit(EstValue.Rho_real(:,1), brainmask_cut_test));
EstMapReal.Rho2_map = flip(img_from_fit(EstValue.Rho_real(:,2), brainmask_cut_test));
EstMapImag.Rho1_map = flip(img_from_fit(EstValue.Rho_imag(:,1), brainmask_cut_test));
EstMapImag.Rho2_map = flip(img_from_fit(EstValue.Rho_imag(:,2), brainmask_cut_test));

figure('Position',[100 100 1200 300]);
subplot(141)
imagesc(EstMapReal.Rho1_map,[0 20]); axis off; axis tight; axis equal;
subplot(142)
imagesc(EstMapReal.Rho2_map,[0 10]); axis off; axis tight; axis equal;
subplot(143)
imagesc(EstMapImag.Rho1_map,[-20 20]); axis off; axis tight; axis equal;
subplot(144)
imagesc(EstMapImag.Rho2_map,[-10 10]); axis off; axis tight; axis equal;
