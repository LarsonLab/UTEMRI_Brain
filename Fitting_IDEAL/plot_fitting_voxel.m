function [] = plot_fitting_voxel(S0, TE, S_fit, iter, rmse_Iter, rmse_NL, Nacq, reslim)

if nargin < 7
    Nacq = 1;
end

if nargin < 8
    reslim = [-0.02 0.02];
end

Necho = length(TE)/Nacq;

set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',14);

figure('Position',[100 100 400 600]);
subplot(311)
for ii = 1:Nacq
    TEind = (ii-1)*Necho+1:ii*Necho;
    set(gca,'ColorOrderIndex',1)
    plot(TE(TEind),real(S0(TEind)),'*'); hold on;
    plot(TE(TEind),imag(S0(TEind)),'o');
    set(gca,'ColorOrderIndex',1)
    plot(TE(TEind),real(S_fit(TEind)),'-');
    plot(TE(TEind),imag(S_fit(TEind)),'--');
end
legend({'Real','Imaginary'},'box','off','location','best');
xlabel('Echo time (s)');
title('Acquired signal');

subplot(312)
for ii = 1:Nacq
    TEind = (ii-1)*Necho+1:ii*Necho;
    set(gca,'ColorOrderIndex',1)
    plot(TE(TEind),real(S0(TEind))-real(S_fit(TEind)),'-'); hold on;
    plot(TE(TEind),imag(S0(TEind))-imag(S_fit(TEind)),':');
end
legend({'Real','Imaginary'},'box','off','location','best');
xlabel('Echo time (s)');
ylim(reslim);
title('IDEAL fitted signal residual');

subplot(313)
plot(1:iter,rmse_Iter); hold on;
plot([1 iter], [rmse_NL rmse_NL], ':');
legend({['IDEAL ' num2str(rmse_Iter(end),4)],['Nonlinear ' num2str(rmse_NL,4)]},...
    'box','off','location','best');
xlabel('Iteration number');
title('RSME of iterations');

end