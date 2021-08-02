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