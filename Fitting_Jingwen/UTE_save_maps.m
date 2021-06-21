function fit_maps = UTE_save_maps(fitting_result, model_name, brainmask, I, AIC_thresh)

modelfit = cat(1,fitting_result.(model_name));
AIC = [modelfit.AIC];

if nargin < 5 || isempty(AIC_thresh)
    datamask = 1;
else
    datamask = (AIC < AIC_thresh);
end

Ncomp = size(modelfit,2);
fit_maps = cell(1,Ncomp);

% initialize maps
for ii = 1:2 % assume 2 comp
    % rho
    fit_maps{ii}.rho = zeros(size(brainmask));
    fit_maps{ii}.rho_err = zeros(size(brainmask));
    fit_maps{ii}.rho_sd = zeros(size(brainmask));
    fit_maps{ii}.rho_cv = zeros(size(brainmask));
    fit_maps{ii}.rho_frac = zeros(size(brainmask));
    % T2
    fit_maps{ii}.T2 = zeros(size(brainmask));
    fit_maps{ii}.T2_err = zeros(size(brainmask));
    fit_maps{ii}.T2_sd = zeros(size(brainmask));
    fit_maps{ii}.T2_cv = zeros(size(brainmask));
    % T1
    fit_maps{ii}.T1 = zeros(size(brainmask));
    fit_maps{ii}.T1_err = zeros(size(brainmask));
    fit_maps{ii}.T1_sd = zeros(size(brainmask));
    fit_maps{ii}.T1_cv = zeros(size(brainmask));
    % df
    fit_maps{ii}.df = zeros(size(brainmask));
    fit_maps{ii}.df_err = zeros(size(brainmask));
    % phi
    for nn = 1:length(modelfit(1,1).phi)
        str = sprintf('phi%i',nn);
        str_err = sprintf('phi%i_err',nn);
        fit_maps{ii}.(str) = zeros(size(brainmask));
        fit_maps{ii}.(str_err) = zeros(size(brainmask));
    end
end

for ii = 1:Ncomp
    % rho
    fit_maps{ii}.rho(I) = mean(cat(1,modelfit(:,ii).rho),2);
    fit_maps{ii}.rho_err(I) = mean(cat(1,modelfit(:,ii).rho_ub) - cat(1,modelfit(:,ii).rho_lb),2);
    fit_maps{ii}.rho_sd = fit_maps{ii}.rho_err/2/1.96;
    fit_maps{ii}.rho_cv = fit_maps{ii}.rho_sd./fit_maps{ii}.rho*100;
    % T2
    fit_maps{ii}.T2(I) = [modelfit(:,ii).T2];
    fit_maps{ii}.T2_err(I) = [modelfit(:,ii).T2_ub] - [modelfit(:,ii).T2_lb];
    fit_maps{ii}.T2_sd = fit_maps{ii}.T2_err/2/1.96;
    fit_maps{ii}.T2_cv = fit_maps{ii}.T2_sd./fit_maps{ii}.T2*100;
    % T1
    if ~isempty(strfind(model_name,'T1'))
        fit_maps{ii}.T1(I) = [modelfit(:,ii).T1];
        fit_maps{ii}.T1_err(I) = [modelfit(:,ii).T1_ub] - [modelfit(:,ii).T1_lb];
        fit_maps{ii}.T1_sd = fit_maps{ii}.T1_err/2/1.96;
        fit_maps{ii}.T1_cv = fit_maps{ii}.T1_sd./fit_maps{ii}.T1*100;
    end
    % phi
    phiArray = cat(1,modelfit(:,ii).phi);
    phiErrArray = cat(2,modelfit(:,ii).phi_ub) - cat(2,modelfit(:,ii).phi_lb);
    for nn = 1:length(modelfit(1,1).phi)
        str = sprintf('phi%i',nn);
        str_err = sprintf('phi%i_err',nn);
        fit_maps{ii}.(str)(I) = phiArray(:,nn);
        fit_maps{ii}.(str_err)(I) = phiErrArray(nn,:);
    end
    
    % SD = 95%CI/2/1.96
    % CV = SD/Estimate*100
    
end

if Ncomp == 2
    % rho fracction
    fit_maps{2}.rho_frac = fit_maps{2}.rho./fit_maps{1}.rho;
    % df
    fit_maps{2}.df(I) = [modelfit(:,2).df];
    fit_maps{2}.df_err(I) = [modelfit(:,2).df_ub] - [modelfit(:,2).df_lb];
end

% AIC
fit_maps{1}.AIC = zeros(size(brainmask));
fit_maps{1}.AIC(I) = AIC;
    
end