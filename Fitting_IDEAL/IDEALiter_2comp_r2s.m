function [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
    IDEALiter_2comp_r2s(S0, TE, phi_iter, R2s_iter, chemshift)

Necho = length(TE);
Ncomp = 2;

% demodulating B0 field
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];

% estmate rho given B0
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE);
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE);
C(:,2) = cos(2*pi*chemshift*TE).*exp(-abs(R2s_iter(2))*TE);
D(:,2) = sin(2*pi*chemshift*TE).*exp(-abs(R2s_iter(2))*TE);

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

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
rhoI_delta = y(2:2:end);

phi_iter = phi_iter + phi_delta;
phi_final = phi_iter;
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end) = rhoR_final;
rho_final(2:2:end) = rhoI_final;

% demodulating B0 field again
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];

% calculate residual: subtract estimated signal from acquired signal
Sreal = diag(repmat([rhoR_final;rhoI_final]',[Necho,1])*[C';-D']);
Simag = diag(repmat([rhoR_final;rhoI_final]',[Necho,1])*[D';C']);
Ssubt = Sdemod - [Sreal; Simag];

% estimate delta R2s
E = zeros(Necho*2, Ncomp);
for ee = 1:Necho
    for nn = 1:Ncomp
        E(ee,nn) = -rhoR_final(nn)*C(ee,nn)*TE(ee) + rhoI_final(nn)*D(ee,nn)*TE(ee);
        E(Necho+ee,nn) = -rhoR_final(nn)*D(ee,nn)*TE(ee) - rhoI_final(nn)*C(ee,nn)*TE(ee);
    end
end

z = (E'*E)\E'*Ssubt;

R2s_delta = z;
R2s_iter = abs(R2s_iter + R2s_delta);
R2s_final = R2s_iter;

% calculate fitted signal and residual
C(:,1) = cos(2*pi*0*TE).*exp(-R2s_final(1)*TE);
D(:,1) = sin(2*pi*0*TE).*exp(-R2s_final(1)*TE);
C(:,2) = cos(2*pi*chemshift*TE).*exp(-R2s_final(2)*TE);
D(:,2) = sin(2*pi*chemshift*TE).*exp(-R2s_final(2)*TE);

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

Sdemod_final = A*rho_final;
S_fit = Sdemod_final(1:Necho) + 1i*Sdemod_final(Necho+1:2*Necho);
S_fit = S_fit.*exp(1i*2*pi*phi_final*TE);

end