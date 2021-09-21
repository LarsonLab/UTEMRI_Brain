function [fitparam, time_IDEAL, rmse_Iter, S_fit, iter] = ...
    IDEALfunc(S0all, TEall, FAall, TR, multipeak, initCond, iterMax, thres)

if nargin < 8
    thres = 0.01;
end

if nargin < 7
    iterMax = 1000;
end

phi_delta = Inf;
R2s_delta = Inf;
R1_delta = Inf;

phi_iter = initCond.phi;
R2s_iter = initCond.R2s;
R1_iter = initCond.R1;

tic;

iter = 0;
rmse_Iter = [];
while (abs(phi_delta) > thres || max(abs(R2s_delta)) > thres || max(abs(R1_delta)) > thres) && iter < iterMax
    
    iter = iter + 1;
    
    [~, ~, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
        IDEALiter_2comp_r2sr1_v3(S0all, TEall, TR, FAall, phi_iter, R2s_iter, R1_iter, multipeak);
    
    phi_iter = phi_final;
    R2s_iter = R2s_final;
    R1_iter = R1_final;
    
    S_res = [real(S0all) - real(S_fit); imag(S0all) - imag(S_fit)];
    
    rmse_Iter(iter) = rms(S_res);
    
end

% fprintf('Iter %i \n', iter);

time_IDEAL = toc;

fitparam.phi = phi_final;
fitparam.R2s = R2s_final;
fitparam.R1 = R1_final;
fitparam.rho = rho_final;
