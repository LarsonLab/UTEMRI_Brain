function [fit_result1, AIC1, fit_result2, AIC2, fit_result2m, AIC2m, fit_result3, AIC3, I] = ute_brain_multiscan_fitting(TE, imall, B0, phi_RF, savefname)

if nargin < 3 || isempty(B0)
    B0 = 3;
end

if nargin < 4 || isempty(phi_RF)
    phi_RF = 0;
end

if nargin < 5 || isempty(savefname)
    writedata_flag = 0;
else
    writedata_flag = 1;
end


im_ute = imall(:,:,:,1);

imsize = size(im_ute);
Nechoes = size(imall, 4);

% only fit non-zero (can be used for masking) values
[I] = find(im_ute);
NI = length(I);

Snormall = max(im_ute(:));

%%
imfit = reshape(imall, [prod(imsize) Nechoes]) / Snormall;
imfit = imfit(I,:);

TEin = TE*1e-3; % ms

parfor Ix = 1:NI
    
    Sin = imfit(Ix,:);
    [fit_result1(Ix,:), AIC1(Ix), fit_result2(Ix,:), AIC2(Ix), fit_result2m(Ix,:), AIC2m(Ix), fit_result3(Ix,:), AIC3(Ix)] = ...
        utebrain_multiscan_fitting_function(TEin, Sin, B0, phi_RF);
    
    %    disp([int2str(Ix) ' of ' int2str(NI)])
end

% repeat with nearby median values to try and eliminate bad fits?

if writedata_flag
    save(savefname,'fit_result1','AIC1', 'fit_result2', 'AIC2', 'fit_result2m', 'AIC2m', 'fit_result3', 'AIC3', 'I', 'imsize');
end




