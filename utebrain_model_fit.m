function [fit_result, rmse, AIC, TEfit, Sfit] = utebrain_model_fit(TE,S,fit_params, general_opts);

default_opts = struct('plot_flag', 0, 'B0', 7, 'num_components', 2, 'complex_fit', 0); % global parameters

names = fieldnames(default_opts);
for k = 1:length(names)
    if ~isfield(general_opts, names(k))
        general_opts.(names{k}) = default_opts.(names{k});
    end
end


% scaling
Snorm = max(abs(S));
S = S/Snorm;

ppm_freq = general_opts.B0*42.57e-3; % kHz
ub_df_default = 4.5*ppm_freq;  % frequency is reversed
lb_df_default = -2*ppm_freq;

X0 = []; lb = []; ub = [];
for n = 1:general_opts.num_components
    X0 = [X0, fit_params(n).rho.est, fit_params(n).T2.est, fit_params(n).df.est, fit_params(n).phi.est];
    
    if isfield(fit_params(n).T2, 'lb')
        lb_T2 = fit_params(n).T2.lb;
    elseif n == 1
        lb_T2 = TE(2) - TE(1);
    else
        lb_T2 = (fit_params(n).T2.est + fit_params(n-1).T2.est)/2;  % or geometric mean
    end
    if isfield(fit_params(n).df, 'lb')
        lb_df = fit_params(n).df.lb;
    else
        lb_df = lb_df_default;
    end
    if isfield(fit_params(n).phi, 'lb')
        lb_phi = fit_params(n).phi.lb;
    else
        lb_phi = fit_params(n).phi.est-2*pi;
    end
    %
    lb = [lb, 0, lb_T2,  lb_df, lb_phi];
    
    if isfield(fit_params(n).T2, 'ub')
        ub_T2 = fit_params(n).T2.ub;
    elseif n == general_opts.num_components
        ub_T2 = 10*TE(end);
    else
        ub_T2 = (fit_params(n).T2.est + fit_params(n+1).T2.est)/2;  % or geometric mean
    end
    if isfield(fit_params(n).df, 'ub')
        ub_df = fit_params(n).df.ub;
    else
        ub_df = ub_df_default;
    end
    if isfield(fit_params(n).phi, 'ub')
        ub_phi = fit_params(n).phi.ub;
    else
        ub_phi = fit_params(n).phi.est+2*pi;
    end
    %
    ub = [ub, 2, ub_T2,  ub_df, ub_phi];
    
end

    function res = model_diff(x)
        
        Sest = utebrain_signal_model(x, general_opts.num_components, TE);
        if general_opts.complex_fit
            res = [real(Sest(:) - S(:)); imag(Sest(:) - S(:))];
        else
            res = abs(Sest(:)) - abs(S(:));
        end
        
    end

%opts = optimset('MaxIter', 2000, 'MaxFunEvals', 2000, 'TolFun', 10^(-7), 'TolX', 10^(-7));
lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);

[X,resnorm,residual,exitflag] = lsqnonlin(@model_diff, X0, lb, ub, lsq_opts);
%exitflag

fit_result = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{});
% T2s = X(4*[0:general_opts.num_components-1] + 2);
% [~,IT2s] = sort(T2s);
for n = 1:general_opts.num_components
%    In = 4*(IT2s(n)-1);
    In = 4*(n-1);
    fit_result(n).rho = X(In + 1)*Snorm;
    fit_result(n).T2 = X(In + 2);
    fit_result(n).df = X(In + 3);
    fit_result(n).phi = X(In + 4);
end

rmse =sqrt(resnorm);

numParams = length(X);
% if general_opts.complex_fit == 0
%     numParams = numParams - 2;  % only relative phases & frequencies can be fit by model
% end
AIC = aic(residual, length(S), numParams);


TEfit=linspace(0,max(TE));
Sfit = utebrain_signal_model(X, general_opts.num_components, TEfit)*Snorm;

if general_opts.plot_flag==1
    Sfit_TE = utebrain_signal_model(X, general_opts.num_components, TE);
    figure
    if general_opts.complex_fit
%     plot(TE,abs(S),'+',TEfit,abs(Sfit/Snorm))
%     subplot(212)
     subplot(311)
    plot(TE,real(S),'b+',TEfit,real(Sfit/Snorm), 'b--')
    subplot(312)
        plot(TE,imag(S),'g+',TEfit,imag(Sfit/Snorm), 'g--')
     subplot(313)
     Sresidual = S(:) - Sfit_TE(:);
    plot(TE,real(Sresidual),'b', TE, imag(Sresidual), 'g')

    else
     subplot(211)
      plot(TE,abs(S),'+',TEfit,abs(Sfit/Snorm))
      subplot(212)
    plot(TE,abs(S(:)) - abs(Sfit_TE(:)))
     
    end
    
end

end