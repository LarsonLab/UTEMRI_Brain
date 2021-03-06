function [fit_result1, AIC1, fit_result2, AIC2, fit_result2m, AIC2m, fit_result3, AIC3] = ...
    utebrain_t1_fitting_function(TEin_all, Sin_all, flips, TR, B0, phi_est, plot_flag)

if nargin < 7
    plot_flag = 0;
end

num_scans = length(Sin_all);

general_opts.plot_flag = plot_flag;  general_opts.B0 = B0;
methylene_freq_est = 3.5* B0*42.57e-3; % kHz

% 1 component complex

general_opts.num_components = 1; general_opts.complex_fit = 1;

fit_params = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});
fit_params_t1 = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});

fit_params(1).rho.est = 1*ones(1,num_scans);
fit_params(1).T2.est = 15;
fit_params(1).df.est = 0;
fit_params(1).phi.est = 0*ones(1,num_scans);

% fit long T2 component
general_opts.use_weights = 0;
general_opts.plot_flag = 0;

[fit_result1, rmse1, AIC1, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_all,fit_params, general_opts);
general_opts.plot_flag = plot_flag;

general_opts.use_weights = 0;

% remove long-T2 phase and frequency (easier to see...)
for n= 1:num_scans
    Sin_corrected{n} = Sin_all{n}(:) .* exp(1i* (2*pi*fit_result1(1).df .* TEin_all{n}(:) - fit_result1(1).phi(n)) );
end

% assume first component fit is good
fit_params(1).T2.est = fit_result1(1).T2;
fit_params(1).T2.lb = fit_result1(1).T2 - 15;
fit_params(1).T2.ub = fit_result1(1).T2 + 15;
fit_params(1).df.est = 0;
fit_params(1).df.lb = -.05;
fit_params(1).df.ub = .05;
fit_params(1).phi.est = 0*ones(1,num_scans);
fit_params(1).phi.lb = -0.05;
fit_params(1).phi.ub = 0.05;


    fit_params1 = fit_params;
    general_opts1 = general_opts;

% fit long T2 component
% [fit_result1_2, rmse1_2, AIC1_2, TEfit, Sfit] = utebrain_model_fit(TEin,Sin_corrected,fit_params, general_opts);

if 1
    % Option 1: add 2 component - add ultrashort, then 3rd component
    
    % add second component
    general_opts.num_components = 2;
    
    fit_params(2).rho.est = 0.1*ones(1,num_scans);
    fit_params(2).T2.est = .5;
    fit_params(2).T2.lb = .1;fit_params(2).T2.ub = 50;
    fit_params(2).df.est = methylene_freq_est;  % strong influence...
    fit_params(2).phi.est = phi_est*ones(1,num_scans);  % dphi from RF pulse
    %fit_params(2).phi.lb = phi_RF(3)-dphi_bound; fit_params(2).phi.ub = phi_RF(3)+dphi_bound;
    
    %   IuT2 = 1:length(Sin_all);
    
    if 1
        % perform magnitude fit as well
        general_opts.complex_fit = 0;
        general_opts.plot_flag = 0;
        [fit_result2m, rmse2m, AIC2m, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_corrected,fit_params, general_opts);
        for n = 1:general_opts.num_components
            fit_params(n).rho.est = fit_result2m(n).rho;
            fit_params(n).T2.est = fit_result2m(n).T2;
            fit_params(n).df.est = fit_result2m(n).df;
        end
        
        general_opts.complex_fit = 1;
        general_opts.plot_flag = plot_flag;
    end
    
    general_opts.plot_flag = 0;
    [fit_result2, rmse2, AIC2, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_corrected,fit_params, general_opts);
    
    % remove long-T2 phase and frequency again
    for n= 1:num_scans
        Sin_corrected{n} = Sin_corrected{n}(:) .* exp(1i* (2*pi*fit_result2(1).df .* TEin_all{n}(:) - fit_result2(1).phi(n)) );
    end
    
    general_opts.plot_flag = plot_flag;
    if plot_flag
        general_opts.use_weights = 0;
        [fit_result1, rmse1, AIC1, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_corrected,fit_params1, general_opts1);
        
        general_opts.use_weights = 0;
    end
    
    for n = 1
        fit_params_t1(n).rho.est = mean(fit_result1(n).rho)* sin(mean(flips));
        fit_params_t1(n).T2.est = fit_result1(n).T2;
        fit_params_t1(n).T2.lb = min(.95* fit_result1(n).T2, 1.05*fit_result1(n).T2);
        fit_params_t1(n).T2.ub = max(.95* fit_result1(n).T2, 1.05*fit_result1(n).T2);
        fit_params_t1(n).df.est = 0;
        fit_params_t1(n).df.lb = -.05;
        fit_params_t1(n).df.ub = +.05;
        fit_params_t1(n).phi.est = 0*ones(1,num_scans);
        fit_params_t1(n).phi.lb = - .05*ones(1,num_scans);
        fit_params_t1(n).phi.ub = + .05*ones(1,num_scans);
    end
    
    % fit absolute fractional component value to each flip angle (3/19)
    % use stripped down version
    % sanity check in data to make sure curves look expected
    
    fit_params_t1(1).T1.est = 0.8;
    
    general_opts1.plot_flag = 0;
    
    [fit_result1, rmse1, AIC1, TEfit, Sfit] = utebrain_t1_model_fit(TEin_all,Sin_corrected, flips, TR, fit_params_t1, general_opts1);
    general_opts1.use_weights = 0;
    
    fit_params_t1(1).T1.est = fit_result1.T1;
    
    [fit_result2, rmse2, AIC2, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_corrected,fit_params, general_opts);
    
    for n = 1:general_opts.num_components
        fit_params(n).rho.est = fit_result2(n).rho;
        fit_params(n).T2.est = fit_result2(n).T2;
        fit_params(n).df.est = fit_result2(n).df;
        fit_params(n).phi.est = fit_result2(n).phi;
    end
    
    
    
    for n = 1:general_opts.num_components
        fit_params_t1(n).rho.est = mean(fit_result2(n).rho)* sin(mean(flips));
        fit_params_t1(n).T2.est = fit_result2(n).T2;
        fit_params_t1(n).T2.lb = min(.95* fit_result2(n).T2, 1.05*fit_result2(n).T2);
        fit_params_t1(n).T2.ub = max(.95* fit_result2(n).T2, 1.05*fit_result2(n).T2);
        fit_params_t1(n).df.est = fit_result2(n).df;
        fit_params_t1(n).df.lb = fit_result2(n).df-.05;
        fit_params_t1(n).df.ub = fit_result2(n).df+.05;
        fit_params_t1(n).phi.est = fit_result2(n).phi;
        fit_params_t1(n).phi.lb = fit_result2(n).phi - .05;
        fit_params_t1(n).phi.ub = fit_result2(n).phi+ .05;
    end
    
  
    
    %fit_params_t1(1).T1.est = .8;
    fit_params_t1(2).T1.est = .3;
    
    
    general_opts.use_weights = 0;
    
    [fit_result2, rmse2, AIC2, TEfit, Sfit] = utebrain_t1_model_fit(TEin_all,Sin_corrected, flips, TR, fit_params_t1, general_opts);
    
    general_opts.use_weights = 0;
    
    general_opts.num_components = 3;
    
    fit_params(3).rho.est = 0.1*ones(1,num_scans);
    fit_params(3).T2.est = 8;
    fit_params(3).T2.lb = .10;fit_params(3).T2.ub = 50;
    fit_params(3).df.est = 0;
    fit_params(3).phi.est = 0*ones(1,num_scans);
    
    [fit_result3, rmse3, AIC3, TEfit, Sfit] = utebrain_multiscan_model_fit(TEin_all,Sin_corrected,fit_params, general_opts);
    
    for n = 1:general_opts.num_components
        fit_params_t1(n).rho.est = mean(fit_result3(n).rho)* sin(mean(flips));
        fit_params_t1(n).T2.est = fit_result3(n).T2;
        fit_params_t1(n).T2.lb = min(.95* fit_result3(n).T2, 1.05*fit_result3(n).T2);
        fit_params_t1(n).T2.ub = max(.95* fit_result3(n).T2, 1.05*fit_result3(n).T2);
        fit_params_t1(n).df.est = fit_result3(n).df;
        fit_params_t1(n).df.lb = fit_result3(n).df-.05;
        fit_params_t1(n).df.ub = fit_result3(n).df+.05;
        fit_params_t1(n).phi.est = fit_result3(n).phi;
        fit_params_t1(n).phi.lb = fit_result3(n).phi-.05;
        fit_params_t1(n).phi.ub = fit_result3(n).phi+.05;
    end
    fit_params_t1(1).T1.est = fit_result2(1).T1; fit_params_t1(2).T1.est = fit_result2(2).T1;  fit_params_t1(3).T1.est = .6;
    
    [fit_result3, rmse3, AIC3, TEfit, Sfit] = utebrain_t1_model_fit(TEin_all,Sin_corrected, flips, TR, fit_params_t1, general_opts);
    
else
    % Option 2: add med component, then ultrashort
    % 2 component - addmed component
    general_opts.num_components = 2;
    
    fit_params(2).rho.est = 0.1;
    fit_params(2).T2.est = 8;
    fit_params(2).T2.lb = .1;fit_params(2).T2.ub = 50;
    fit_params(2).df.est = methylene_freq_est;
    fit_params(2).phi.est = 0;
    
    [fit_result2, rmse2, AIC2, TEfit, Sfit] = utebrain_model_fit(TEin,Sin_corrected,fit_params, general_opts);
    
    % second component is now good
    for n = 1:general_opts.num_components
        fit_params(n).rho.est = fit_result2(n).rho;
        fit_params(n).T2.est = fit_result2(n).T2;
        fit_params(n).df.est = fit_result2(n).df;
        fit_params(n).phi.est = fit_result2(n).phi;
    end
    % tight on bounds too?
    
    general_opts.num_components = 3;
    
    fit_params(3).rho.est = 0.1;
    fit_params(3).T2.est = .5;
    fit_params(3).T2.lb = .100;fit_params(3).T2.ub = 5;
    fit_params(3).df.est = methylene_freq_est;  % strong influence...
    fit_params(3).phi.est = phi_est;  % dphi from RF pulse
    
    [fit_result3, rmse3, AIC3, TEfit, Sfit] = utebrain_model_fit(TEin,Sin_corrected,fit_params, general_opts);
    
end
