function [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
    IDEALiter_r2s(S0, TE, phi_iter, R2s_iter, A, C, D)

Necho = length(TE);

% demodulating B0 field and R2* decay
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE + R2s_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];

% estmate rho given B0
rho = (A'*A)\A'*Sdemod;
rhoR = rho(1:2:end);
rhoI = rho(2:2:end);

% calculate residual: subtract estimated signal from acquired signal
Sreal = diag(repmat([rhoR;rhoI]',[Necho,1])*[C';-D']);
Simag = diag(repmat([rhoR;rhoI]',[Necho,1])*[D';C']);
Ssubt = Sdemod - [Sreal; Simag];

% estimate delta phi and delta rho
GR = 2*pi*TE.*diag(repmat([rhoR;rhoI]',[Necho,1])*[-D';-C']);
GI = 2*pi*TE.*diag(repmat([rhoR;rhoI]',[Necho,1])*[C';-D']);
RR = TE.*diag(repmat([rhoR;rhoI]',[Necho,1])*[-C';D']);
RI = TE.*diag(repmat([rhoR;rhoI]',[Necho,1])*[-D';-C']);

B = [[RR;RI] [GR;GI] A];

y = (B'*B)\B'*Ssubt;

R2s_delta = y(1);
phi_delta = y(2);
rhoR_delta = y(3:2:end);
rhoI_delta = y(4:2:end);

phi_iter = phi_iter + phi_delta;
R2s_iter = R2s_iter + R2s_delta;

% calculate fitted signal and residual
phi_final = phi_iter;
R2s_final = R2s_iter;
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end) = rhoR_final;
rho_final(2:2:end) = rhoI_final;

Sdemod_final = A*rho_final;
S_fit = Sdemod_final(1:Necho) + 1i*Sdemod_final(Necho+1:2*Necho);
S_fit = S_fit.*exp(1i*2*pi*phi_final*TE - R2s_final*TE);

end