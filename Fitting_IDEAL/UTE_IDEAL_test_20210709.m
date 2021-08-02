clc; clear;

%% add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils'));
addpath(genpath('/home/jyao3/020_UTE_Brain'));

%% load data
load('test_IDEAL_data.mat');

figure('Position',[100 100 400 400]);
for ii = 1:3
    h(ii) = plot(TEin{ii},real(Sin{ii}),'-'); hold on;
    plot(TEin{ii},imag(Sin{ii}),':');
end
legend(h, {'FA = 6deg','12deg','18deg'},'box','off','location','best');
xlabel('Echo time (ms)');
% export_fig(sprintf('Svoxel'),'-png','-transparent'); close;

%% set up parameters

B0 = 3; % Tesla
TR = 9.1e-3; % s
methylene_freq_est = 3.5 * B0 * 42.57e-3; % kHz

%% NL fit for initialization

tic;

general_opts.B0 = B0;
general_opts.use_weights = 0;
general_opts.TR = TR;
general_opts.num_components = 1;
general_opts.complex_fit = 1;

fit_params = struct('rho',{}, 'T2',{}, 'df', {}, 'phi',{}, 'T1', {});
fit_params(1).T2.est = 15;
fit_params(1).T1.est = 1.5;
fit_params(1).T2.lb = .1;
fit_params(1).T2.ub = 50;
fit_params(1).df.est = 0;
fit_params(1).phi.est = 0;
fit_params(1).rho.est = 1/sin(mean(aflips));

[fitting_result, rmse, AIC, TEfit, Sfit, Sfit_TE] = ...
    UTE_T1T2_model_fit_singlePhi(TEin, Sin, aflips, TR, fit_params, general_opts);
fitting_result.RMSE = rmse;
fitting_result.AIC = AIC;

% plot_fitting(general_opts, TEin, Sin, TEfit, Sfit, Sfit_TE);
% subplot(211); title('1 comp T2+T1');
% export_fig(sprintf('S_nl1'),'-png','-transparent'); close;

% remove long-T2 phase and frequency
num_scans = length(TEin);
for n= 1:num_scans
    Sin_corrected{n} = Sin{n}(:) .* ...
        exp(1i* (2*pi*fitting_result(1).df .* TEin{n}(:)) );
end

general_opts.num_components = 2;
general_opts.plot_flag = 0;
general_opts.fixPhi = 0;
general_opts.fixDf = 1;
general_opts.phi_RF = 0;
general_opts.methylene_freq_est = methylene_freq_est;
[fitting_result_step2, TEfit, Sfit, Sfit_TE] = ...
    UTE_fitting_function_step2(TEin, Sin, aflips, general_opts, fitting_result, [0 0 0]);

% plot_fitting(general_opts, TEin, Sin_corrected, TEfit, Sfit, Sfit_TE);
% subplot(211); title('2 comp T2+T1');
% export_fig(sprintf('S_nl2'),'-png','-transparent'); close;

rmse_NL = rms([real(Sin_corrected{1})-real(Sfit_TE{1}); ...
    imag(Sin_corrected{1})-imag(Sfit_TE{1})]);

time_NL = toc;

%% print results

fprintf('== Nonlinear fitting == \n');
fprintf('rho 1 %f, rho 2 %f \n', ...
    fitting_result_step2.comp2T1(1).rho, fitting_result_step2.comp2T1(2).rho);
fprintf('rho fraction %f \n', ...
    fitting_result_step2.comp2T1(2).rho/fitting_result_step2.comp2T1(1).rho);
fprintf('T2s 1 %f ms, T2s 2 %f ms \n', ...
    fitting_result_step2.comp2T1(1).T2, fitting_result_step2.comp2T1(2).T2);
fprintf('T1 1 %f s, T1 2 %f s \n', ...
    fitting_result_step2.comp2T1(1).T1, fitting_result_step2.comp2T1(2).T1);
fprintf('B0 %f Hz \n', ...
    fitting_result.df + fitting_result_step2.comp2T1(1).df);

%% setting up IDEAL

phi_init = 60; % fitting_result.df; % Hz
R2s_init = 30; % 1000/fitting_result.T2; % s-1

chemshift = -methylene_freq_est; % 3.5 * B0 * 42.57e-3

S0 = transpose(Sin{1});
TE = transpose(TEin{1})/1000; % s

S0all = transpose([Sin{:}]);
TEall = transpose([TEin{:}])/1000; % s

Necho = length(TE);
Ncomp = 2;

%% original IDEAL

C(:,1) = cos(2*pi*0*TE);
D(:,1) = sin(2*pi*0*TE);
C(:,2) = cos(2*pi*chemshift*TE);
D(:,2) = sin(2*pi*chemshift*TE);

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

phi_delta = Inf;
phi_iter = phi_init;

iter = 0;
rmse_Iter = [];

while abs(phi_delta) > 0.001

    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    [Sdemod, Ssubt, S_fit, phi_final, rho_final, phi_delta] = ...
        IDEALiter_ori(S0, TE, phi_iter, A, C, D);
    
    phi_iter = phi_iter + phi_delta;
    
    S_res = [real(S0) - S_fit(1:Necho); imag(S0) - S_fit(Necho+1:2*Necho)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

plot_fitting_voxel(S0, TE, S_fit, iter, rmse_Iter, rmse_NL);

% export_fig(sprintf('IDEALori'),'-png','-transparent'); close;

%% IDEAL with R2*

C(:,1) = cos(2*pi*0*TE);
D(:,1) = sin(2*pi*0*TE);
C(:,2) = cos(2*pi*chemshift*TE);
D(:,2) = sin(2*pi*chemshift*TE);

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

phi_delta = Inf;
R2s_delta = Inf;
phi_iter = phi_init;
R2s_iter = R2s_init;

iter = 0;
rmse_Iter = [];

while abs(phi_delta) > 0.001 || abs(R2s_delta) > 0.001
    
    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    
    [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
        IDEALiter_r2s(S0, TE, phi_iter, R2s_iter, A, C, D);
    
    phi_iter = phi_iter + phi_delta;
    R2s_iter = R2s_iter + R2s_delta;
    
    S_res = [real(S0) - real(S_fit); imag(S0) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

plot_fitting_voxel(S0, TE, S_fit, iter, rmse_Iter, rmse_NL);

% export_fig(sprintf('IDEAL_R2s'),'-png','-transparent'); close;

%% IDEAL with multi-compartment R2*

tic;

phi_delta = Inf;
R2s_delta = Inf;
phi_iter = 63;
R2s_iter = [27; 2987];

iter = 0;
rmse_Iter = [];

while abs(phi_delta) > 0.01 || max(abs(R2s_delta)) > 0.01
    
    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    
    [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
        IDEALiter_2comp_r2s(S0, TE, phi_iter, R2s_iter, chemshift);
    
    phi_iter = phi_iter + phi_delta;
    R2s_iter = R2s_iter + R2s_delta;
    
    S_res = [real(S0) - real(S_fit); imag(S0) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

time_IDEAL = toc;

plot_fitting_voxel(S0, TE, S_fit, iter, rmse_Iter, rmse_NL);

% export_fig(sprintf('IDEAL_R2s_2comp'),'-png','-transparent'); close;

rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
rho_fract = abs(rho_mag(2))/(sum(abs(rho_mag)));

fprintf('== IDEAL fitting 2comp == \n');
fprintf('rho 1 %f, rho 2 %f \n', abs(rho_mag));
fprintf('rho fraction %f \n', rho_fract);
fprintf('T2s 1 %f ms, T2s 2 %f ms \n', 1000./R2s_final);
fprintf('B0 %f Hz \n', phi_final);

%% IDEAL with multi-compartment R2* with multi-peak

tic;

phi_delta = Inf;
R2s_delta = Inf;
phi_iter = 65;
R2s_iter = [30; 3000];

% from 10.1073/PNAS.1115107109
multipeak.alpha = [0.743 0.124 0.111 0.021];
multipeak.alpha = multipeak.alpha/sum(multipeak.alpha);
multipeak.chemshift = [1.55 0.9 3.2 1.3]-4.8; % ppm
multipeak.chemshift = multipeak.chemshift * B0 * 42.57e-3;

iter = 0;
rmse_Iter = [];

S0 = S0all(1:8);
TE = TEall(1:8);

while abs(phi_delta) > 0.01 || max(abs(R2s_delta)) > 0.01
    
    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    
    [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
        IDEALiter_2comp_r2s_multipeak(S0, TE, phi_iter, R2s_iter, multipeak);
    
    phi_iter = phi_iter + phi_delta;
    R2s_iter = R2s_iter + R2s_delta;
    
    S_res = [real(S0) - real(S_fit); imag(S0) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

time_IDEAL = toc;

plot_fitting_voxel(S0, TE, S_fit, iter, rmse_Iter, rmse_NL);

% pause(0.1); export_fig(sprintf('IDEAL_R2s_2comp_multipeak'),'-png','-transparent'); close;

rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
rho_fract = abs(rho_mag(2))/(sum(abs(rho_mag)));

fprintf('== IDEAL fitting 2comp == \n');
fprintf('rho 1 %f, rho 2 %f \n', abs(rho_mag));
fprintf('rho fraction %f \n', rho_fract);
fprintf('T2s 1 %f ms, T2s 2 %f ms \n', 1000./R2s_final);
fprintf('B0 %f Hz \n', phi_final);

%% IDEAL with multi-compartment R2*, multi-peak, multi FA

tic;

phi_delta = Inf;
R2s_delta = Inf;
phi_iter = 100;
R2s_iter = [1; 1000];

% from 10.1073/PNAS.1115107109
multipeak.alpha = [0.743 0.124 0.111 0.021];
multipeak.alpha = multipeak.alpha/sum(multipeak.alpha);
multipeak.chemshift = [1.55 0.9 3.2 1.3]-4.8; % ppm
multipeak.chemshift = multipeak.chemshift * B0 * 42.57e-3;

iter = 0;
rmse_Iter = [];

S0test = S0all(1:24);
TEtest = TEall(1:24);

while (abs(phi_delta) > 0.01 || max(abs(R2s_delta)) > 0.01) && iter < 1000
    
    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    
    [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
        IDEALiter_2comp_r2s_multipeak(S0test, TEtest, phi_iter, R2s_iter, multipeak);
    
    phi_iter = phi_iter + phi_delta;
    R2s_iter = R2s_iter + R2s_delta;
    
    S_res = [real(S0test) - real(S_fit); imag(S0test) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

time_IDEAL = toc;

Nacq = 3;
plot_fitting_voxel(S0test, TEtest, S_fit, iter, rmse_Iter, rmse_NL, Nacq);

% pause(0.1); export_fig(sprintf('IDEAL_R2s_2comp_multiFA'),'-png','-transparent'); close;

rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
rho_fract(1) = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));
rho_fract(2) = abs(rho_mag(4))/(sum(abs(rho_mag(3:4))));
rho_fract(3) = abs(rho_mag(6))/(sum(abs(rho_mag(5:6))));

fprintf('== IDEAL fitting 2comp == \n');
fprintf('rho 1 %f, rho 2 %f \n', abs(rho_mag));
fprintf('rho fraction %f \n', rho_fract);
fprintf('T2s 1 %f ms, T2s 2 %f ms \n', 1000./R2s_final);
fprintf('B0 %f Hz \n', phi_final);

%% calculate T1 from rho

TR = 9.1e-3; % s

rho_mag = abs(rho_final(1:2:end) + 1i*rho_final(2:2:end));
rho_mag_comp1 = rho_mag(1:2:end);
rho_mag_comp2 = rho_mag(2:2:end);

ind = [1,2;2,3;1,3];
T1_comp1 = zeros(1,3);
T1_comp2 = zeros(1,3);
for ii = 1:size(ind,1)
    FA1 = aflips(ind(ii,1)); FA2 = aflips(ind(ii,2));
    rho1 = rho_mag_comp1(ind(ii,1)); rho2 = rho_mag_comp1(ind(ii,2));
    E1 = (sin(FA1)*rho2 - sin(FA2)*rho1)/(cos(FA2)*sin(FA1)*rho2 - cos(FA1)*sin(FA2)*rho1);
    T1_comp1(ii) = -TR/log(E1); % s
    rho1 = rho_mag_comp2(ind(ii,1)); rho2 = rho_mag_comp2(ind(ii,2));
    E1 = (sin(FA1)*rho2 - sin(FA2)*rho1)/(cos(FA2)*sin(FA1)*rho2 - cos(FA1)*sin(FA2)*rho1);
    T1_comp2(ii) = -TR/log(E1); % s
end

fprintf('comp 1 T1 (s) %f %f %f mean %f \n', T1_comp1, mean(T1_comp1));
fprintf('comp 2 T1 (s) %f %f %f mean %f \n', T1_comp2, mean(T1_comp2));

%% IDEAL with multi-compartment R2*+R1, multi-peak

tic;

phi_delta = Inf;
R2s_delta = Inf;
phi_iter = 60;
R2s_iter = [1; 1000];
R1_iter = [0.5; 3];

% from 10.1073/PNAS.1115107109
multipeak.alpha = [0.743 0.124 0.111 0.021];
multipeak.alpha = multipeak.alpha/sum(multipeak.alpha);
multipeak.chemshift = [1.55 0.9 3.2 1.3]-4.8; % ppm
multipeak.chemshift = multipeak.chemshift * B0 * 42.57e-3;

iter = 0;
rmse_Iter = [];

FAall = repmat(aflips,[Necho,1]); FAall = FAall(:);

while (abs(phi_delta) > 0.001 || max(abs(R2s_delta)) > 0.001 || max(abs(R1_delta)) > 0.001) && iter < 100
    
    iter = iter + 1;
    fprintf('Iter %i \n', iter);
    
    [Sdemod, Ssubt, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
        IDEALiter_2comp_r2sr1_multipeak(S0all, TEall, TR, FAall, phi_iter, R2s_iter, R1_iter, multipeak);
    
    phi_iter = phi_iter + phi_delta;
    R2s_iter = R2s_iter + R2s_delta;
    R1_iter = R1_iter + R1_delta;
    
    S_res = [real(S0all) - real(S_fit); imag(S0all) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

time_IDEAL = toc;

Nacq = 3;
plot_fitting_voxel(S0all, TEall, S_fit, iter, rmse_Iter, rmse_NL, Nacq);

% pause(0.1); export_fig(sprintf('IDEAL_R2sR1_2comp_multiFA'),'-png','-transparent'); close;

rho_mag = rho_final(1:2:end) + 1i*rho_final(2:2:end);
rho_fract = abs(rho_mag(2))/(sum(abs(rho_mag(1:2))));

fprintf('== IDEAL fitting 2comp == \n');
fprintf('rho 1 %f, rho 2 %f \n', abs(rho_mag));
fprintf('rho fraction %f \n', rho_fract);
fprintf('T2s 1 %f ms, T2s 2 %f ms \n', 1000./R2s_final);
fprintf('T1 1 %f s, T1 2 %f s \n', 1./R1_final);
fprintf('B0 %f Hz \n', phi_final);
