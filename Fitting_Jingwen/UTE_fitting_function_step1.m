function [fitting_result] = ...
    UTE_fitting_function_step1(TEin_all, Sin_all, flips, TR, B0, phi_est, plot_flag)

if nargin < 7
    plot_flag = 0;
end

num_scans = length(Sin_all);

general_opts.plot_flag = plot_flag;
general_opts.B0 = B0;
general_opts.use_weights = 0;

Sin_corrected = Sin_all;

% 1 component complex

general_opts.num_components = 1;
general_opts.complex_fit = 1;

fit_params = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});
fit_params(1).rho.est = 1*ones(1,num_scans);
fit_params(1).T2.est = 15;
fit_params(1).T1.est = 1.5;
fit_params(1).T2.lb = .1;
fit_params(1).T2.ub = 50;
fit_params(1).df.est = 0;
fit_params(1).phi.est = 0*phi_est*ones(1,num_scans);
fit_params(1).rho.est = 1/sin(mean(flips));

% 1 component + T1

[fitting_result.comp1T1, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
    UTE_T1T2_model_fit(TEin_all, Sin_corrected, flips, TR, fit_params, general_opts);
fitting_result.comp1T1.RMSE = rmse;
fitting_result.comp1T1.AIC = AIC;

if plot_flag
    plot_fitting(general_opts, TEin_all, Sin_corrected, TEfit, Sfit, Sfit_TE);
    subplot(211); title('1 comp T2+T1');
    export_fig('fitplot_Comp1wT1','-png','-transparent'); close;
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
        subplot(211)
        plot(TE_all{n},abs(S_all{n}/Snorm),'+','Color', ColorOrder(n,:)); hold on;
        l(n) = plot(TEfit{n},abs(Sfit{n}/Snorm),'--','Color', ColorOrder(n,:));
        xlabel('TE (ms)'); ylabel('Normalized signal');
        ylim([-1 1]);
        subplot(212)
        plot(TE_all{n},abs(S_all{n}(:)/Snorm) - abs(Sfit_TE{n}(:)/Snorm),'-','Color', ColorOrder(n,:));
        xlabel('TE (ms)'); ylabel('Residual signal');
        ylim([-0.1 0.1]);
    end
end
subplot(211)
legend(l,{'FA = 6deg','12deg','18deg'},'box','off','location','best');

end