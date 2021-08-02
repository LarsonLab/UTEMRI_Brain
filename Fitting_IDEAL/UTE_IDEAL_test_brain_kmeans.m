clc; clear;

%% add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/home/jyao3/020_UTE_Brain'));

%% load data
load('test_IDEAL_slice.mat');

figure('Position',[100 100 400 400]);
imagesc(flip(abs(im_norm_test{1}(:,:,1,1)))); axis off; axis tight; axis equal;

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
    parfor ii = 1:indN
        Smean = transpose(mean(Sinput(idx == ii,:),1));
        FAmean = transpose(mean(FAinput(idx == ii,:),1));
        
        if m == 1
            phi_iter = 0; % Hz
            R2s_iter = [30; 3000]; % s-1
            R1_iter = [0.5; 3]; % s-1
        else
            phi_iter = nanmean(EstValue.Phi(idx == ii),1);
            R2s_iter = nanmean(EstValue.R2s(idx == ii,:),1)';
            R1_iter = nanmean(EstValue.R1(idx == ii,:),1)';
        end
        
        % IDEAL iterations
        phi_delta = Inf;
        R2s_delta = Inf;
        R1_delta = Inf;
        iter = 0;
        
        while (abs(phi_delta) > 1 || max(abs(R2s_delta)) > 0.1 || max(abs(R1_delta)) > 0.1) && iter < 100
            iter = iter + 1;
            %         fprintf('Iter %i \n', iter);
            [~, ~, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
                IDEALiter_2comp_r2sr1_multipeak(Smean, TEall, TR, FAmean, phi_iter, R2s_iter, R1_iter, multipeak);
            phi_iter = phi_final;
            R2s_iter = R2s_final;
            R1_iter = R1_final;
        end
        
        % assign component according to R2
        [R2s_iter,indComp] = sort(R2s_iter);
        R1_iter = R1_iter(indComp);
        
        % save estimation values into arrays
        rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
        Rho(ii,:) = abs(rho_mag)';
        Rho_real(ii,:) = real(rho_mag)';
        Rho_imag(ii,:) = imag(rho_mag)';
        R2s(ii,:) = R2s_iter';
        R1(ii,:) = R1_iter';
        Phi(ii) = phi_iter;
        Rho_frac(ii) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));
        
        S_res = [real(Smean) - real(S_fit); imag(Smean) - imag(S_fit)];
        Rmse(ii) = rms(S_res);
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
    
end

% if mod(Ix,100) == 0
%     disp(['Processing voxel ' num2str(Ix) ' of ' num2str(length(I)) ' (' num2str(Ix/length(I)*100,2) '%)']);
% end

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
