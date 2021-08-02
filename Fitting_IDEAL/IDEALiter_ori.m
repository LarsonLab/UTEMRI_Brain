function [Sdemod, Ssubt, S_fit, phi_final, rho_final, phi_delta] = ...
    IDEALiter_ori(S0, TE, phi_iter, A, C, D)

Necho = length(TE);

% demodulating B0 field and R2* decay
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
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

B = [[GR;GI] A];

y = (B'*B)\B'*Ssubt;

phi_delta = y(1);
rhoR_delta = y(2:2:end);
rhoI_delta = y(3:2:end);

phi_iter = phi_iter + phi_delta;

% calculate fitted signal and residual
phi_final = phi_iter;
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end) = rhoR_final;
rho_final(2:2:end) = rhoI_final;

Sdemod_final = A*rho_final;
S_fit = Sdemod_final(1:Necho) + 1i*Sdemod_final(Necho+1:2*Necho);
S_fit = S_fit.*exp(1i*2*pi*phi_final*TE);
S_fit = [real(S_fit); imag(S_fit)];

end