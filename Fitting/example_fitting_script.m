clear variables;
writedata_flag = 0;

%%
Ifolder = 1;
foldername{1} = '../SampleData/';  % must have git lfs ...
foldername{1} = '/data/larson/brain_uT2/2018-07-20-3T-DTI-volunteer/';
B0 = 3;

%% load data
%load([foldername{Ifolder} '/ute_32echo_random-csreconallec_l2_r0p01.mat'])
load([foldername{Ifolder} '/P26112.7-csreconallec_l2_r0p01.mat'])

%% parameter estimate
methylene_freq_est = 3.5* B0*42.57e-3; % kHz
phi_RF = .0562;

%% whole volume
load([foldername{Ifolder} '/brainmask.mat'])

[I] = find(brainmask);
NI = length(I);
imsize = size(imall(:,:,:,1));
Nechoes = size(imall, 4);

fitsize = [1 NI];

im_masked = imall(:,:,:,1) .* brainmask;
Snormall = max(im_masked(:));

%%
for plusfit = 0:1;
    
    
    if plusfit
        imfit = reshape(imallplus, [prod(imsize) Nechoes]) / Snormall;
    else
        imfit = reshape(imall, [prod(imsize) Nechoes]) / Snormall;
    end
    imfit = imfit(I,:);
    
    TEin = TE*1e-3; % ms
    
    parfor Ix = 1:NI
        
        Sin = imfit(Ix,:);
        [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
            utebrain_fitting_function(TEin, Sin, B0, phi_RF);
        
        %    disp([int2str(Ix) ' of ' int2str(NI)])
    end
    
    % repeat with nearby median values to try and eliminate bad fits?
    
    if plusfit
        savefname = 'fit_results_plus';
    else
        savefname = 'fit_results';
    end
    if writedata_flag
        save(sprintf('%s/%s', foldername{Ifolder}, savefname),'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
    end
end



