clc; clear;

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

sliceN = 50;

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

% Sinput = []; % Nvoxel x 24
% for i = 1:3
%     Sinput = [Sinput imfit{i}];
% end
% 
% B1_scale  = FAmapfit/AFI_nominal_flip;
% FAall = repmat(flips,[Necho/length(TEin),1]); FAall = FAall(:);
% FAinput = repmat(FAall',[length(I),1]) .* repmat(B1_scale,[1,Necho]); % Nvoxel x 24

%% IDEAL

Rho = zeros(length(I),2);
Rho_real = zeros(length(I),2);
Rho_imag = zeros(length(I),2);
R2s = zeros(length(I),2);
R1 = zeros(length(I),2);
Phi = zeros(length(I),1);
Rho_frac = zeros(length(I),1);
Rmse = zeros(length(I),1);
Iter = zeros(length(I),1);

initCond = cell(1,length(I));
parfor Ix = 1:100 % 1:length(I)
    
    if mod(Ix,100) == 0
        disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
    end
    
    initCond{Ix}.phi = 60; % Hz
    initCond{Ix}.R2s = [30; 3000]; % s-1
    initCond{Ix}.R1 = [0.5; 3]; % s-1
    
    S0all = transpose(Sinput(Ix,:));
    FAall = transpose(FAinput(Ix,:));

%     Sin = cell(1,3);
%     for i = 1:3
%         Sin{i} = imfit{i}(Ix,:);
%     end
%     S0all = transpose([Sin{:}]);
%     
%     B1_scale  = FAmapfit(Ix) ./ AFI_nominal_flip;
%     aflips = flips .* B1_scale;
%     FAall = repmat(aflips,[Necho/length(TEin),1]); FAall = FAall(:);
    
    % IDEAL iterations
    iterMax = 200;
    thres = 0.1;
    [fitparam, ~, rmse_Iter, S_fit, iter] = ...
        IDEALfunc(S0all, TEall, FAall, TR, multipeak, initCond{Ix}, iterMax, thres);
    
    if isnan(fitparam.phi)
        initCond{Ix}.phi = 60;
        initCond{Ix}.R2s = [30; 3000]; % s-1
        initCond{Ix}.R1 = [0.5; 3]; % s-1
        [fitparam, ~, rmse_Iter, S_fit, iter] = ...
            IDEALfunc(S0all, TEall, FAall, TR, multipeak, initCond{Ix}, iterMax, thres);
    end
       
    % assign component according to R2
    [fitparam.R2s,indComp] = sort(fitparam.R2s);
    fitparam.R1 = fitparam.R1(indComp);
    
    % save estimation values into arrays
    rho_mag = fitparam.rho(1:2:end) + 1i*fitparam.rho(2:2:end);
    rho_mag = rho_mag(indComp);
    
    Rho(Ix,:) = abs(rho_mag);
    Rho_real(Ix,:) = real(rho_mag);
    Rho_imag(Ix,:) = imag(rho_mag);
    R2s(Ix,:) = fitparam.R2s;
    R1(Ix,:) = fitparam.R1;
    Phi(Ix) = fitparam.phi;
    Rho_frac(Ix) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));
    
    Rmse(Ix) = rmse_Iter(end);
    Iter(Ix) = iter;
    
end

% plot_fitting_voxel(S0all, TEall, S_fit, iter, rmse_Iter, 0, 3)

%% plot maps

Rho1_map = flip(img_from_fit(Rho(:,1), mask_slice>0));
Rho2_map = flip(img_from_fit(Rho(:,2), mask_slice>0));
RhoF_map = flip(img_from_fit(Rho_frac(:,1), mask_slice>0));
T2s1_map = flip(img_from_fit(1000./R2s(:,1), mask_slice>0));
T2s2_map = flip(img_from_fit(1000./R2s(:,2), mask_slice>0));
T11_map = abs(flip(img_from_fit(1./R1(:,1), mask_slice>0)));
T12_map = abs(flip(img_from_fit(1./R1(:,2), mask_slice>0)));
Phi_map = flip(img_from_fit(Phi(:,1), mask_slice>0));
Rmse_map = flip(img_from_fit(Rmse(:,1), mask_slice>0));
Iter_map = flip(img_from_fit(Iter(:,1), mask_slice>0));

figure('Position',[100 100 1600 600]);
subplot(241)
imagesc(Rho1_map,[0 20]); axis off; axis tight; axis equal; title('Rho1'); colorbar;
subplot(242)
imagesc(Rho2_map,[0 10]); axis off; axis tight; axis equal; title('Rho2'); colorbar;
subplot(243)
imagesc(RhoF_map,[0 0.3]); axis off; axis tight; axis equal; title('Rho fraction'); colorbar;
subplot(244)
imagesc(Phi_map,[-100 200]); axis off; axis tight; axis equal; title('delta B0'); colorbar;

subplot(247)
imagesc(T2s1_map,[0 100]); axis off; axis tight; axis equal; title('T2s comp1'); colorbar;
subplot(248)
imagesc(T2s2_map,[0 1]); axis off; axis tight; axis equal; title('T2s comp2'); colorbar;
subplot(245)
imagesc(T11_map,[0 5]); axis off; axis tight; axis equal; title('T1 comp1'); colorbar;
subplot(246)
imagesc(T12_map,[0 1]); axis off; axis tight; axis equal; title('T1 comp2'); colorbar;
% export_fig(sprintf('IDEAL_slice'),'-png','-transparent'); close;

figure('Position',[100 100 1200 300]);
subplot(141)
imagesc(Rmse_map,[0 0.02]); axis off; axis tight; axis equal; title('RMSE'); colorbar;
subplot(142)
imagesc(Iter_map,[0 200]); axis off; axis tight; axis equal; title('iter number'); colorbar;
% export_fig(sprintf('IDEAL_slice2'),'-png','-transparent'); close;

%% plot more

EstMapReal.Rho1_map = flip(img_from_fit(Rho_real(:,1), mask_slice>0));
EstMapReal.Rho2_map = flip(img_from_fit(Rho_real(:,2), mask_slice>0));
EstMapImag.Rho1_map = flip(img_from_fit(Rho_imag(:,1), mask_slice>0));
EstMapImag.Rho2_map = flip(img_from_fit(Rho_imag(:,2), mask_slice>0));

figure('Position',[100 100 1200 300]);
subplot(141)
imagesc(EstMapReal.Rho1_map,[0 20]); axis off; axis tight; axis equal; colorbar; title('rho real comp1'); 
subplot(142)
imagesc(EstMapReal.Rho2_map,[0 10]); axis off; axis tight; axis equal; colorbar; title('rho real comp2'); 
subplot(143)
imagesc(EstMapImag.Rho1_map,[-20 20]); axis off; axis tight; axis equal; colorbar; title('rho imag comp1'); 
subplot(144) 
imagesc(EstMapImag.Rho2_map,[-10 10]); axis off; axis tight; axis equal; colorbar; title('rho imag comp2'); 
% export_fig(sprintf('IDEAL_slice_complex'),'-png','-transparent');
