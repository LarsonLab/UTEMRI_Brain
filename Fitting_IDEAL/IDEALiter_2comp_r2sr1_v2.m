function [Sdemod, Ssubt, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
    IDEALiter_2comp_r2sr1_v2(S0, TE, TR, FA, phi_iter, R2s_iter, R1_iter, multipeak)

Necho = length(TE);

Ncomp = 2;

% demodulating B0 field
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];

% estmate rho given B0, R1, R2*
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

rho = (A'*A)\A'*Sdemod;
rhoR = rho(1:2:end);
rhoI = rho(2:2:end);

% calculate residual: subtract estimated signal from acquired signal
Sfit = A*rho;
Ssubt = Sdemod - Sfit;

% estimate delta phi, delta rho, delta R2*
G1 = 2*pi*repmat(TE,[2 1]).*([-D, -C; C, -D]*[rhoR;rhoI]);
B = G1;

G3 = zeros(size(G1));
for nn = 1:Ncomp
    G3(:,nn) = -repmat(TE,[2 1]).*([C(:,nn), -D(:,nn); D(:,nn), C(:,nn)]*[rhoR(nn);rhoI(nn)]);
    B = [B [C(:,nn) -D(:,nn); D(:,nn) C(:,nn)] G3(:,nn)];
end

y = (B'*B)\B'*Ssubt;

phi_delta = y(1);
phi_iter = phi_iter + phi_delta;
phi_final = phi_iter;

rhoR_delta = y(2:3:end);
rhoI_delta = y(3:3:end);
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end,:) = rhoR_final;
rho_final(2:2:end,:) = rhoI_final; % rho1R rho1I rho2R ...

R2s_delta = y(4:3:end);
R2s_iter = abs(R2s_iter + R2s_delta);
R2s_final = R2s_iter;

% estimate delta r1
Sdemod_C = S0.*exp(-1i*2*pi*phi_final*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];
Sfit = A*rho_final;
Ssubt = Sdemod - Sfit;

C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));

B = zeros(size(G1));
for nn = 1:Ncomp
    E1 = exp(-TR*abs(R1_iter(nn)));
    B(:,nn) = TR*E1/(1-E1)*([C(:,nn), -D(:,nn); D(:,nn), C(:,nn)]*[rhoR_final(nn);rhoI_final(nn)]);
end

z = (B'*B)\B'*Ssubt;

R1_delta = z;
R1_iter = abs(R1_iter + R1_delta);
R1_final = R1_iter;

% calculate fitted data and residual
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(2))));

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

Sdemod_final = A*rho_final;
S_fit = Sdemod_final(1:Necho) + 1i*Sdemod_final(Necho+1:2*Necho);
S_fit = S_fit.*exp(1i*2*pi*phi_final*TE);

end