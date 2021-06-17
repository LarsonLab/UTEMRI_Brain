function [fit_result, rmse, AIC, TEfit, Sfit] = utebrain_multiscan_model_fit(TE_all,S_all,fit_params, general_opts);

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
Snorm = max(abs(S(:)));

for n = 1:num_scans
S_all{n} = S_all{n}/Snorm;
end

ppm_freq = general_opts.B0*42.57e-3; % kHz
ub_df_default = 4.5*ppm_freq;  % frequency is reversed
lb_df_default = -2*ppm_freq;

X0 = []; lb = []; ub = [];
for n = 1:general_opts.num_components
     X0 = [X0, fit_params(n).rho.est, fit_params(n).T2.est, fit_params(n).df.est, fit_params(n).phi.est];  % T2, df assumed to be the same across experiments
   
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
        lb_phi = fit_params(n).phi.lb*ones(1,num_scans);
    else
        lb_phi = fit_params(n).phi.est-2*pi*ones(1,num_scans);
    end
    %
    lb = [lb, 0*ones(1,num_scans), lb_T2,  lb_df, lb_phi];
    
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
        ub_phi = fit_params(n).phi.ub*ones(1,num_scans);
    else
        ub_phi = fit_params(n).phi.est+2*pi*ones(1,num_scans);
    end
    %
    ub = [ub, 2*ones(1,num_scans), ub_T2,  ub_df, ub_phi];
    
end



    function res = model_diff(x)
        res = [];
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
            
            x_scan = generate_x_scan(x,n);
                    Sest{n} = utebrain_signal_model(x_scan, general_opts.num_components, TE_all{n});
        if general_opts.complex_fit
            if general_opts.use_weights == 1;
                %res = [res; weight{n}(:) .* real(Sest{n}(:) - S_all{n}(:)); weight{n}(:) .* imag(Sest{n}(:) - S_all{n}(:))];
            else
                res = [res; real(Sest{n}(:) - S_all{n}(:)); imag(Sest{n}(:) - S_all{n}(:))];
            end
      
            % weight{n}(:) * real(Sest... * imag(Sest...
        else
            res = [res; abs(Sest{n}(:)) - abs(S_all{n}(:))];
            % weight{n}(:) * abs(Sest...
        end
        end
    end

    function x_out = generate_x_scan(x_in, Iscan)
        x_out = [];
        for m = 1:general_opts.num_components
                % x = [(component 1) rho_scan1 rho_scan2 ... T2 df phi_scan1 phi_scan2 ..., (component 2) rho_scan1 rho_scan2 ... T2 df phi_scan1 phi_scan2 , ...]
                    % x_scan = [(component 1) rho_scanm T2 df phi_scanm, (component 2) rho_scanm T2 df phi_scanm, ...]
                    component_offset = (2*num_scans+2)*(m-1);
                x_out = [x_out, x_in(Iscan + component_offset),  x_in(num_scans+1 + component_offset ) ...
                    x_in(num_scans+2 + component_offset ) x_in(num_scans+2+Iscan + component_offset)];
        end
    end

%opts = optimset('MaxIter', 2000, 'MaxFunEvals', 2000, 'TolFun', 10^(-7), 'TolX', 10^(-7));
lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);

[X,resnorm,residual,~,~,~,J] = lsqnonlin(@model_diff, X0, lb, ub, lsq_opts);
  
CI = nlparci(X,residual,'jacobian',J);
X_lb = CI(:,1);
X_ub = CI(:,2);

fit_result = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{});
% T2s = X(4*[0:general_opts.num_components-1] + 2);
% [~,IT2s] = sort(T2s);
for n = 1:general_opts.num_components
    %    In = 4*(IT2s(n)-1);
    component_offset = (2*num_scans+2)*(n-1);
    fit_result(n).rho(1:num_scans) = X(component_offset + [1:num_scans])*Snorm;
    fit_result(n).T2 = X(component_offset + num_scans+1);
    fit_result(n).df = X(component_offset + num_scans+2);
    fit_result(n).phi = X(component_offset + num_scans+2+ [1:num_scans]);
    % lower bound
    fit_result(n).rho_lb(1:num_scans) = X_lb(component_offset + [1:num_scans])*Snorm;
    fit_result(n).T2_lb = X_lb(component_offset + num_scans+1);
    fit_result(n).df_lb = X_lb(component_offset + num_scans+2);
    fit_result(n).phi_lb = X_lb(component_offset + num_scans+2+ [1:num_scans]);
    % upper bound
    fit_result(n).rho_ub(1:num_scans) = X_ub(component_offset + [1:num_scans])*Snorm;
    fit_result(n).T2_ub = X_ub(component_offset + num_scans+1);
    fit_result(n).df_ub = X_ub(component_offset + num_scans+2);
    fit_result(n).phi_ub = X_ub(component_offset + num_scans+2+ [1:num_scans]);
end

rmse =sqrt(resnorm);

numParams = length(X);
% if general_opts.complex_fit == 0
%     numParams = numParams - 2;  % only relative phases & frequencies can be fit by model
% end
AIC = aic(residual, length(S), numParams);

for n =1:num_scans
TEfit{n}=linspace(0,max(TE_all{n}));
            X_scan = generate_x_scan(X,n);

Sfit{n} = utebrain_signal_model(X_scan, general_opts.num_components, TEfit{n})*Snorm;

if general_opts.plot_flag==1
    Sfit_TE = utebrain_signal_model(X_scan, general_opts.num_components, TE_all{n});
%    figure
%     if general_opts.complex_fit
%        %     plot(TE,abs(S),'+',TEfit,abs(Sfit/Snorm))
%        %     subplot(212)
%         subplot(411)
%         plot(TE_all{n},real(S_all{n}),'b+',TEfit{n},real(Sfit{n}/Snorm), 'b--')
%         subplot(412)
%         plot(TE_all{n},real(S_all{n}),'b+',TEfit{n},real(Sfit{n}/Snorm), 'b--')
%         subplot(413)
%         plot(TE_all{n},imag(S_all{n}),'g+',TEfit{n},imag(Sfit{n}/Snorm), 'g--')
%         subplot(414)
%         Sresidual = S_all{n}(:) - Sfit_TE(:);
%         plot(TE_all{n},real(Sresidual),'b', TE_all{n}, imag(Sresidual), 'g')
%         
%     else
%         subplot(211)
%         plot(TE_all{n},abs(S_all{n}),'+',TEfit{n},abs(Sfit{n}/Snorm))
%         subplot(212)
%         plot(TE_all{n},abs(S_all{n}(:)) - abs(Sfit_TE(:)))
%         
%     end
    
end

end
end