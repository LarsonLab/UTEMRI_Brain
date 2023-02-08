function [fitparam, time_NONLIN, rmse, S_fit] = ...
    LSQNONLINfunc(S0all, Acq, multipeak, initCond, lb, ub, tRF, iterMax, Wall)

% Input:
%       S0all: array of complex acquired signal, size Nacquisition x 1
%       Acq.TE: array of TEs (second), size Nacquisition x 1
%       Acq.FA: array of flip angles (rad), B1 corrected, size Nacquisition x 1
%       Acq.TR: array of TRs (second), size Nacquisition x 1
%       multipeak.alpha: multipeak amplitude fractions (normalized to 1) of myelin lipid peaks
%       multipeak.chemshift: multipeak chemical shifts (Hz) of myelin lipid peaks
%   Fitting parameters 
%       (in the order of PDwr: water proton equilibrium signal real
%                        PMmr: myelin proton equilibrium signal real
%                        PDwi: water proton equilibrium signal imaginary
%                        PMmi: myelin proton equilibrium signal imaginary
%                        R2w: water R2* (s-1)
%                        R2m: myelin R2* (s-1)
%                        R1w: water R1 (s-1)
%                        R1m: myelin R1 v
%                        DFw: water frequency shift from B0 (field offset, Hz)
%       initCond: initial condition
%       lb: lower bound
%       ub: upper bound

if nargin < 9
    Wall = ones(size(Acq.TE(:)));
    Wall = [Wall; Wall];
end

if nargin < 8
    iterMax = 500;
end

if nargin < 7
    tRF = 0; % 0.8*10^-3; % s
end

if nargin < 5
        % PDwr PDmr PDwi PDmi R2w  R2m   R1w R1m  DFw
    lb = [-inf -inf -inf -inf 0.5  1000  0.2 0.2  -inf];
    ub = [+inf +inf +inf +inf 100  10000 10  10   +inf];
end

if nargin < 4
              % PDwr PDmr PDwi PDmi R2w R2m  R1w R1m DFw
    initCond = [10   1    0    0    25  1600 0.5 2   30];
end

Cond.init = initCond;
Cond.lb = lb; % [-inf -inf -inf -inf 0.5  1000  0.2 2 -inf];
Cond.ub = ub; % [+inf +inf +inf +inf 100  10000 10  2 +inf];

tic;

lsq_opts = optimset('Display','none','MaxIter',iterMax,'MaxFunEvals',iterMax); % 'iter'
[X,~,~,~,~,~,~] = lsqnonlin(@cost_comp1, Cond.init, Cond.lb, Cond.ub, lsq_opts);

time_NONLIN = toc;

fitparam.phi = X(9);
fitparam.R2s = [X(5) X(6)];
fitparam.R1 = [X(7) X(8)];
fitparam.rho = [X(1)+1i*X(3) X(2)+1i*X(4)];

[S_res, S_fit] = cost_comp1(X);
rmse = rms(S_res);

%% cost function
    function [res, S] = cost_comp1(param)
        
        PDw = param(1)+1i*param(3);
        PDm = param(2)+1i*param(4);
        R2w = param(5);
        R2m = param(6);
        R1w = param(7);
        R1m = param(8);
        DFw = param(9);
        
        tRF = 0.8*10^-3; % s
        tau_w = R2w*tRF/2;
        tau_m = R2m*tRF/2;
        
        S = zeros(length(Acq.TE(:)),1);
        for nn = 1:length(Acq.TE(:))
            phi_w = sqrt(Acq.FA(nn)^2-tau_w^2);
            phi_m = sqrt(Acq.FA(nn)^2-tau_m^2);
            
            S(nn) = PDw .* ... % M0
                sin(phi_w) .* (1-exp(-Acq.TR(nn)*R1w)) .* exp(-tau_w) .* (Acq.FA(nn)/phi_w) ...
                ./ (1-exp(-Acq.TR(nn)*R1w) .* exp(-tau_w) .* (cos(phi_w) + tau_w/phi_w*sin(phi_w))) .*  ...  %T1 weighting
                exp(-Acq.TE(nn)*R2w) .* ... % T2 star weighting
                exp(1i*2*pi*DFw*Acq.TE(nn));  % frequency shift;
            for jj = 1:length(multipeak.alpha)
                S(nn) = S(nn) + ...
                    PDm*multipeak.alpha(jj) .* ... % M0
                    sin(phi_m) .* (1-exp(-Acq.TR(nn)*R1m)) .* exp(-tau_m) .* (Acq.FA(nn)/phi_m) ...
                    ./ (1-exp(-Acq.TR(nn)*R1m) .* exp(-tau_m) .* (cos(phi_m) + tau_m/phi_m*sin(phi_m))) .*  ...  %T1 weighting
                    exp(-Acq.TE(nn)*R2m) .* ... % T2 star weighting
                    exp(1i*2*pi*(DFw + multipeak.chemshift(jj))*Acq.TE(nn)); % chemical shift + frequency shift
            end
            
%             S(nn) = PDw .* ... % M0
%                 sin(Acq.FA(nn)) .* (1-exp(-Acq.TR(nn)*R1w)) ./ (1-cos(Acq.FA(nn))*exp(-Acq.TR(nn)*R1w)) .*  ...  %T1 weighting
%                 exp(-Acq.TE(nn)*R2w) .* ... % T2 star weighting
%                 exp(1i*2*pi*DFw*Acq.TE(nn));  % frequency shift;
%             for jj = 1:length(multipeak.alpha)
%                 S(nn) = S(nn) + ...
%                     PDm*multipeak.alpha(jj) .* ... % M0
%                     sin(Acq.FA(nn)) .* (1-exp(-Acq.TR(nn)*R1m)) ./ (1-cos(Acq.FA(nn))*exp(-Acq.TR(nn)*R1m)) .*  ...  %T1 weighting
%                     exp(-Acq.TE(nn)*R2m) .* ... % T2 star weighting
%                     exp(1i*2*pi*(DFw + multipeak.chemshift(jj))*Acq.TE(nn)); % chemical shift + frequency shift
%             end
        end
        
        res = abs([real(S - S0all); imag(S - S0all)].*Wall);
        
    end

end