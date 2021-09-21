clc; clear;

%% add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/home/jyao3/020_UTE_Brain'));

%% load data
load('test_IDEAL_slice.mat');

figure('Position',[100 100 400 400]);
imagesc(flip(abs(im_norm_test{1}(:,:,1,1)))); axis off; axis tight; axis equal;
% export_fig(sprintf('SliceImg'),'-png','-transparent'); close;

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

Sinput = []; % Nvoxel x 24
for i = 1:3
    Sinput = [Sinput imfit{i}];
end

B1_scale  = FAmapfit/AFI_nominal_flip;
FAall = repmat(flips,[Necho/length(TEin),1]); FAall = FAall(:);
FAinput = repmat(FAall',[length(I),1]) .* repmat(B1_scale,[1,Necho]); % Nvoxel x 24

%% IDEAL

flagT1IDEAL = 1;

Rho = zeros(length(I),2);
Rho_real = zeros(length(I),2);
Rho_imag = zeros(length(I),2);
R2s = zeros(length(I),2);
R1 = zeros(length(I),2);
Phi = zeros(length(I),1);
Rho_frac = zeros(length(I),1);
Rmse = zeros(length(I),1);
Iter = zeros(length(I),1);

for Ix = 31 % 1:100 % 1:length(I)
    
    if mod(Ix,100) == 0
        disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
    end
    
    Sin = cell(1,3);
    for i = 1:3
        Sin{i} = imfit{i}(Ix,:);
    end
    S0all = transpose([Sin{:}]);
    
    B1_scale  = FAmapfit(Ix) ./ AFI_nominal_flip;
    aflips = flips .* B1_scale;
    FAall = repmat(aflips,[Necho/length(TEin),1]); FAall = FAall(:);
       
    phi_iter = 60; % Hz
    R2s_iter = [30; 3000]; % s-1
    R1_iter = [0.5; 3]; % s-1
    
%     phi_delta = Inf;
%     R2s_delta = Inf;
%     iter = 0;
%     while (abs(phi_delta) > 0.1 || max(abs(R2s_delta)) > 0.1) && iter < 200
%         iter = iter + 1;
%         [~, ~, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
%                 IDEALiter_2comp_r2s_multipeak(S0all, TEall, phi_iter, R2s_iter, multipeak);
%         phi_iter = phi_final;
%         R2s_iter = R2s_final;
%     end
%     
%     rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
%     rho_mag_comp1 = rho_mag(1:2:end);
%     rho_mag_comp2 = rho_mag(2:2:end);
%     
%     % calculate T1 according to rho
%     ind = [1,2;2,3;1,3];
%     T1_comp1 = zeros(1,3);
%     T1_comp2 = zeros(1,3);
%     for ii = 1:size(ind,1)
%         FA1 = aflips(ind(ii,1)); FA2 = aflips(ind(ii,2));
%         rho1 = abs(rho_mag_comp1(ind(ii,1))); rho2 = abs(rho_mag_comp1(ind(ii,2)));
%         E1 = (sin(FA1)*rho2 - sin(FA2)*rho1)/(cos(FA2)*sin(FA1)*rho2 - cos(FA1)*sin(FA2)*rho1);
%         T1_comp1(ii) = -TR/log(E1); % s
%         rho1 = abs(rho_mag_comp2(ind(ii,1))); rho2 = abs(rho_mag_comp2(ind(ii,2)));
%         E1 = (sin(FA1)*rho2 - sin(FA2)*rho1)/(cos(FA2)*sin(FA1)*rho2 - cos(FA1)*sin(FA2)*rho1);
%         T1_comp2(ii) = -TR/log(E1); % s
%     end
%     R1_iter = abs(1./[mean(T1_comp1(3)); mean(T1_comp2(3))]);
    
    phi_delta = Inf;
    R2s_delta = Inf;
    R1_delta = Inf;
    iter = 0;
    while (abs(phi_delta) > 0.01 || max(abs(R2s_delta)) > 0.1 || max(abs(R1_delta)) > 0.1) && iter < 200
        iter = iter + 1;
        [~, ~, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
            IDEALiter_2comp_r2sr1_v3(S0all, TEall, TR, FAall, phi_iter, R2s_iter, R1_iter, multipeak);
        R1_iter = R1_final;
        phi_iter = phi_final;
        R2s_iter = R2s_final;
    end
    
    rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
    
    % assign component according to R2
    [R2s_iter,indComp] = sort(R2s_iter);
    R1_iter = R1_iter(indComp);
    rho_mag = rho_mag(indComp);
    
    Rho(Ix,:) = abs(rho_mag);
    Rho_real(Ix,:) = real(rho_mag(1:2));
    Rho_imag(Ix,:) = imag(rho_mag(1:2));
    R2s(Ix,:) = R2s_iter;
    R1(Ix,:) = R1_iter;
    Phi(Ix) = phi_iter;
    Rho_frac(Ix) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));
    
    S_res = [real(S0all) - real(S_fit); imag(S0all) - imag(S_fit)];
    Rmse(Ix) = rms(S_res);
    Iter(Ix) = iter;
    
end

plot_fitting_voxel(S0all, TEall, S_fit, iter, rms(S_res), 0, 3)

%% plot maps

Rho1_map = flip(img_from_fit(Rho(:,1), brainmask_cut_test));
Rho2_map = flip(img_from_fit(Rho(:,2), brainmask_cut_test));
RhoF_map = flip(img_from_fit(Rho_frac(:,1), brainmask_cut_test));
T2s1_map = flip(img_from_fit(1000./R2s(:,1), brainmask_cut_test));
T2s2_map = flip(img_from_fit(1000./R2s(:,2), brainmask_cut_test));
T11_map = abs(flip(img_from_fit(1./R1(:,1), brainmask_cut_test)));
T12_map = abs(flip(img_from_fit(1./R1(:,2), brainmask_cut_test)));
Phi_map = flip(img_from_fit(Phi(:,1), brainmask_cut_test));
Rmse_map = flip(img_from_fit(Rmse(:,1), brainmask_cut_test));
Iter_map = flip(img_from_fit(Iter(:,1), brainmask_cut_test));

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

EstMapReal.Rho1_map = flip(img_from_fit(Rho_real(:,1), brainmask_cut_test));
EstMapReal.Rho2_map = flip(img_from_fit(Rho_real(:,2), brainmask_cut_test));
EstMapImag.Rho1_map = flip(img_from_fit(Rho_imag(:,1), brainmask_cut_test));
EstMapImag.Rho2_map = flip(img_from_fit(Rho_imag(:,2), brainmask_cut_test));

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
