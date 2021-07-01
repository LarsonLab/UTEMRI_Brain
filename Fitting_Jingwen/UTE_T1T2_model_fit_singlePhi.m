function [fit_result, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
    UTE_T1T2_model_fit_singlePhi(TE_all,S_all,flips, TR, fit_params, general_opts)

default_opts = struct('plot_flag', 0, 'B0', 3, 'num_components', 2, 'complex_fit', 0); % global parameters

names = fieldnames(default_opts);
for k = 1:length(names)
    if ~isfield(general_opts, names(k))
        general_opts.(names{k}) = default_opts.(names{k});
    end
end

num_scans = length(S_all);

% scaling
S = [];
for n = 1:num_scans
    S = [S;S_all{n}(:)];
end
Snorm = max(abs(S(:))) * max(sin(flips)) * 10;  % help to equilibrate parameter ranges for search

for n = 1:num_scans
    S_all{n} = S_all{n}/Snorm;
end

ppm_freq = general_opts.B0*42.57e-3; % kHz
ub_df_default = 4.5*ppm_freq;  % frequency is reversed
lb_df_default = -2*ppm_freq;
ub_T1_default = 5; % s
lb_T1_default = 0.1; % s

X0 = []; lb = []; ub = [];
for n = 1:general_opts.num_components
    X0 = [X0, fit_params(n).rho.est, fit_params(n).T2.est, fit_params(n).df.est, ...
        fit_params(n).phi.est(1), fit_params(n).T1.est];  % T2, df, T1 assumed to be the same across experiments
    
    if isfield(fit_params(n).T2, 'lb')
        lb_T2 = fit_params(n).T2.lb;
    elseif n == 1
        lb_T2 = TE_all{1}(2) - TE_all{1}(1);
    else
        lb_T2 = (fit_params(n).T2.est + fit_params(n-1).T2.est)/2;  % or geometric mean
    end
    if isfield(fit_params(n).df, 'lb')
        lb_df = fit_params(n).df.lb;
    else
        lb_df = lb_df_default;
    end
    if isfield(fit_params(n).phi, 'lb')
        lb_phi = fit_params(n).phi.lb(1);
    else
        lb_phi = fit_params(n).phi.est(1)-2*pi;
    end
    if isfield(fit_params(n).T1, 'lb')
        lb_T1 = fit_params(n).T1.lb;
    else
        lb_T1 = lb_T1_default;
    end
    %
    lb = [lb, 0, lb_T2,  lb_df, lb_phi, lb_T1];
    
    if isfield(fit_params(n).T2, 'ub')
        ub_T2 = fit_params(n).T2.ub;
    elseif n == general_opts.num_components
        ub_T2 = 10*TE_all{1}(end);
    else
        ub_T2 = (fit_params(n).T2.est + fit_params(n+1).T2.est)/2;  % or geometric mean
    end
    if isfield(fit_params(n).df, 'ub')
        ub_df = fit_params(n).df.ub;
    else
        ub_df = ub_df_default;
    end
    if isfield(fit_params(n).phi, 'ub')
        ub_phi = fit_params(n).phi.ub(1);
    else
        ub_phi = fit_params(n).phi.est(1)+2*pi;
    end
    if isfield(fit_params(n).T1, 'ub')
        ub_T1 = fit_params(n).T1.ub;
    else
        ub_T1 = ub_T1_default;
    end
    %
    ub = [ub, 20, ub_T2,  ub_df, ub_phi, ub_T1];
    
end

    function [res, resJ] = model_diff(x)
        % resJ
        res = [];
        x_scan = generate_x_scan(x);
        for n = 1:num_scans
            
            if general_opts.use_weights == 1
                % density compensation
                for k = 1:length(TE_all{n}) - 1
                    weight{n}(k) = (TE_all{n}(k) + TE_all{n}(k+1)) / 2;
                end
                weight{n}(length(TE_all{n})) = 2 * (TE_all{n}(length(TE_all{n}))- TE_all{n}(length(TE_all{n}) -1));
                
            else
                weight{n} = ones(1,8);
            end
            
            Sest{n} = UTE_T1T2_signal_model(x_scan, TE_all{n}, flips(n), general_opts);
            if general_opts.complex_fit
                if general_opts.use_weights == 1;
                    res = [res; weight{n}(:) .* real(Sest{n}(:) - S_all{n}(:)); weight{n}(:) .* imag(Sest{n}(:) - S_all{n}(:))];
                else
                    res = [res; real(Sest{n}(:) - S_all{n}(:)); imag(Sest{n}(:) - S_all{n}(:))];
                end
            else
                res = [res; abs(Sest{n}(:)) - abs(S_all{n}(:))];
            end
        end
    end

    function x_out = generate_x_scan(x_in)
        x_out = [];
        for m = 1:general_opts.num_components
            % x = [(component 1) rho T2 df phi_scan1 phi_scan2 ...T1 , (component 2) rho T2 df phi_scan1 phi_scan2 ...T1, ...]
            % x_scan = [(component 1) rho T2 df phi_scanm, T1 (component 2) rho T2 df phi_scanm T1, ...]
            component_offset = 5*(m-1);
            x_out = [x_out, x_in(1 + component_offset),  x_in(2 + component_offset ) ...
                x_in(3 + component_offset ) x_in(4 + component_offset) ...
                x_in(5 + component_offset)];
        end
    end

%opts = optimset('MaxIter', 2000, 'MaxFunEvals', 2000, 'TolFun', 10^(-7), 'TolX', 10^(-7));
lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);

[X,resnorm,residual,~,~,~,J] = lsqnonlin(@model_diff, X0, lb, ub, lsq_opts);
%exitflag
  
CI = nlparci(X,residual,'jacobian',J);
X_lb = CI(:,1);
X_ub = CI(:,2);

fit_result = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});
% T2s = X(4*[0:general_opts.num_components-1] + 2);
% [~,IT2s] = sort(T2s);
for n = 1:general_opts.num_components
    component_offset = 5*(n-1);
    fit_result(n).rho = X(component_offset + 1)*Snorm;
    fit_result(n).T2 = X(component_offset + 2);
    fit_result(n).df = X(component_offset + 3);
    fit_result(n).phi = X(component_offset + 4)*ones(1,num_scans);
    fit_result(n).T1 = X(component_offset + 5);
    % lower bound
    fit_result(n).rho_lb = X_lb(component_offset + 1)*Snorm;
    fit_result(n).T2_lb = X_lb(component_offset + 2);
    fit_result(n).df_lb = X_lb(component_offset + 3);
    fit_result(n).phi_lb = X_lb(component_offset + 4)*ones(1,num_scans);
    fit_result(n).T1_lb = X_lb(component_offset + 5);
    % upper bound
    fit_result(n).rho_ub = X_ub(component_offset + 1)*Snorm;
    fit_result(n).T2_ub = X_ub(component_offset + 2);
    fit_result(n).df_ub = X_ub(component_offset + 3);
    fit_result(n).phi_ub = X_ub(component_offset + 4)*ones(1,num_scans);
    fit_result(n).T1_ub = X_ub(component_offset + 5);
end

rmse =sqrt(resnorm);

numParams = length(X);
% if general_opts.complex_fit == 0
%     numParams = numParams - 2;  % only relative phases & frequencies can be fit by model
% end
AIC = aic(residual, length(S), numParams);

for n = 1:num_scans
    TEfit{n}=linspace(0,max(TE_all{n}));
    X_scan = generate_x_scan(X);
    
    Sfit{n} = utebrain_t1_signal_model(X_scan, general_opts.num_components, TEfit{n}, flips(n), TR)*Snorm;
    Sfit_TE{n} = utebrain_t1_signal_model(X_scan, general_opts.num_components, TE_all{n}, flips(n), TR)*Snorm;
    
end

end