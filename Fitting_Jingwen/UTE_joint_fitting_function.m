function [fitting_result] = ...
    UTE_joint_fitting_function(TEin_all, Sin_all, flips, general_opts)

B0 = general_opts.B0;
TR = general_opts.TR;
plot_flag = general_opts.plot_flag;

num_scans = length(Sin_all);

general_opts.use_weights = 0;
methylene_freq_est = 3.5*B0*42.57e-3; % kHz

Sin_corrected = Sin_all;

% 1 component complex

general_opts.num_components = 1;
if ~isfield(general_opts, 'complex_fit')
    general_opts.complex_fit = 1;
end

fit_params = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});
% comp1
fit_params(1).rho.est = 1/sin(mean(flips));
fit_params(1).T2.est = 15;
fit_params(1).T1.est = 1.5;
fit_params(1).T2.lb = .1;
fit_params(1).T2.ub = 50;
fit_params(1).df.est = 0;
fit_params(1).phi.est = 0*ones(1,num_scans);
if general_opts.fixPhi
    fit_params(1).phi.lb = 0*ones(1,num_scans);
    fit_params(1).phi.ub = 0*ones(1,num_scans);
end
% comp2
fit_params(2).rho.est = 0.1*ones(1,num_scans);
fit_params(2).T2.est = .5;
fit_params(2).T2.lb = .1;
fit_params(2).T2.ub = 50;
fit_params(2).T1.est = 0.3;
if general_opts.fixT1comp2
    fit_params(2).T1.est = general_opts.T1median;
    fit_params(2).T1.lb = general_opts.T1median;
    fit_params(2).T1.ub = general_opts.T1median;
end
fit_params(2).df.est = methylene_freq_est;  % strong influence...
fit_params(2).phi.est = general_opts.phi_RF*ones(1,num_scans);  % dphi from RF pulse
if general_opts.fixPhi
    fit_params(2).phi.est = 0*ones(1,num_scans);
    fit_params(2).phi.lb = 0*ones(1,num_scans);
    fit_params(2).phi.ub = 0*ones(1,num_scans);
end

% 1 component + T1
if general_opts.singlePhi
    [fitting_result.comp1T1, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
        UTE_T1T2_model_fit_singlePhi(TEin_all, Sin_corrected, flips, TR, fit_params, general_opts);
else
    [fitting_result.comp1T1, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
        UTE_T1T2_model_fit(TEin_all, Sin_corrected, flips, TR, fit_params, general_opts);
end
fitting_result.comp1T1.RMSE = rmse;
fitting_result.comp1T1.AIC = AIC;

if plot_flag
    plot_fitting(general_opts, TEin_all, Sin_corrected, TEfit, Sfit, Sfit_TE);
%     subplot(211); title('1 comp T2+T1');
%     export_fig('fitplot_Comp1wT1','-png','-transparent'); close;
end

fit_params(1).T1.est = fitting_result.comp1T1.T1;

% remove long-T2 phase and frequency
for n= 1:num_scans
    Sin_corrected{n} = Sin_corrected{n}(:) .* ...
        exp(1i* (2*pi*fitting_result.comp1T1(1).df .* TEin_all{n}(:) - fitting_result.comp1T1(1).phi(n)) );
end

% 2 component

if ~isfield(general_opts, 'complex_fit')
    general_opts.complex_fit = 1;
end
general_opts.num_components = 2;

% 2 component + T1

fit_params(1).rho.est = fitting_result.comp1T1(1).rho;
fit_params(2).rho.est = 0.1 * fitting_result.comp1T1(1).rho;

if general_opts.fixDf
    fit_params(2).df.lb = general_opts.methylene_freq_est;
    fit_params(2).df.ub = general_opts.methylene_freq_est;
end

if general_opts.singlePhi
    [fitting_result.comp2T1, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
        UTE_T1T2_model_fit_singlePhi(TEin_all,Sin_corrected, flips, TR, fit_params, general_opts);
else
    [fitting_result.comp2T1, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
        UTE_T1T2_model_fit(TEin_all,Sin_corrected, flips, TR, fit_params, general_opts);
end
fitting_result.comp2T1(1).RMSE = rmse;
fitting_result.comp2T1(1).AIC = AIC;

if plot_flag
    plot_fitting(general_opts, TEin_all, Sin_corrected, TEfit, Sfit, Sfit_TE);
%     subplot(211); title('2 comp T2+T1');
%     export_fig('fitplot_Comp2wT1','-png','-transparent'); close;
end

end

function [] = plot_fitting(general_opts, TE_all, S_all, TEfit, Sfit, Sfit_TE)

S = [];
num_scans = length(S_all);
for n = 1:num_scans
    S = [S;S_all{n}(:)];
end
Snorm = max(abs(S(:)));

ColorOrder = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];

figure('Position',[100 100 400 500]);
for n = 1:num_scans
    if general_opts.complex_fit
        subplot(211)
        plot(TE_all{n},real(S_all{n}/Snorm),'*','Color', ColorOrder(n,:)); hold on;
        l(n) = plot(TEfit{n},real(Sfit{n}/Snorm),'--','Color', ColorOrder(n,:));
        plot(TE_all{n},imag(S_all{n}/Snorm),'o','Color', ColorOrder(n,:));
        plot(TEfit{n},imag(Sfit{n}/Snorm),'--','Color', ColorOrder(n,:));
        xlabel('TE (ms)'); ylabel('Normalized signal');
        ylim([-1 1]);
        subplot(212)
        Sresidual = S_all{n}(:)/Snorm - Sfit_TE{n}(:)/Snorm;
        plot(TE_all{n},real(Sresidual),'-','Color', ColorOrder(n,:)); hold on;
        plot(TE_all{n}, imag(Sresidual),'-','Color', ColorOrder(n,:)); hold on;
        xlabel('TE (ms)'); ylabel('Residual signal');
        ylim([-0.1 0.1]);
    else
        subplot(221)
        plot(TE_all{n},abs(S_all{n}/Snorm),'+','Color', ColorOrder(n,:)); hold on;
        l(n) = plot(TEfit{n},abs(Sfit{n}/Snorm),'--','Color', ColorOrder(n,:));
        xlabel('TE (ms)'); ylabel('Normalized signal');
        ylim([-1 1]);
        subplot(222)
        plot(TEfit{n},real(Sfit{n}/Snorm),'-','Color', ColorOrder(n,:)); hold on;
        plot(TEfit{n},imag(Sfit{n}/Snorm),':','Color', ColorOrder(n,:));
        xlabel('TE (ms)'); ylabel('Fitted signal');
        ylim([-1 1]);
        subplot(212)
        plot(TE_all{n},abs(S_all{n}(:)/Snorm) - abs(Sfit_TE{n}(:)/Snorm),'-','Color', ColorOrder(n,:)); hold on;
        xlabel('TE (ms)'); ylabel('Residual signal');
        ylim([-0.1 0.1]);
    end
end
subplot(221)
legend(l,{'FA = 6deg','12deg','18deg'},'box','off','location','best');

end