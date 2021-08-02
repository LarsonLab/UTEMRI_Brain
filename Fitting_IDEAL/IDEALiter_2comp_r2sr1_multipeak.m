function [Sdemod, Ssubt, S_fit, phi_final, R2s_final, R1_final, rho_final, phi_delta, R2s_delta, R1_delta] = ...
    IDEALiter_2comp_r2sr1_multipeak(S0, TE, TR, FA, phi_iter, R2s_iter, R1_iter, multipeak)

Necho = length(TE);
Ncomp = 2;

% demodulating B0 field
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)];

% estmate rho given B0
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_iter(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(2))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(2))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));

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

% estimate delta phi and delta rho
G = 2*pi*repmat(TE,[2 1]).*([-D, -C; C, -D]*[rhoR;rhoI]);
B = [G A];

y = (B'*B)\B'*Ssubt;

phi_delta = y(1);
rhoR_delta = y(2:2:end);
rhoI_delta = y(3:2:end);

phi_final = phi_iter + phi_delta;
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end,:) = rhoR_final;
rho_final(2:2:end,:) = rhoI_final;

% demodulating B0 field again and calculate the residual
Sdemod_C = S0.*exp(-1i*2*pi*phi_final*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)]; % Necho x Nflipangle*2 Real
Sfit = A*rho_final;
Ssubt = Sdemod - Sfit;

% estimate delta R2s
E = zeros(Necho*2,Ncomp);
for nn = 1:Ncomp
    E(:,nn) = -repmat(TE,[2 1]).*([C(:,nn), -D(:,nn); D(:,nn), C(:,nn)]*[rhoR_final(nn);rhoI_final(nn)]);
end
% for ee = 1:Necho
%     for nn = 1:Ncomp
%         E(ee,nn) = -rhoR_final(nn)*C(ee,nn)*TE(ee) + rhoI_final(nn)*D(ee,nn)*TE(ee);
%         E(Necho+ee,nn) = -rhoR_final(nn)*D(ee,nn)*TE(ee) - rhoI_final(nn)*C(ee,nn)*TE(ee);
%     end
% end

z = (E'*E)\E'*Ssubt;

R2s_delta = z;
R2s_final = abs(R2s_iter + R2s_delta);

% calculate residual given updated R2*
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(1))))./(1-cos(FA)*exp(-TR*abs(R1_iter(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(2))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_iter(2))))./(1-cos(FA)*exp(-TR*abs(R1_iter(2))));

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

Sfit = A*rho_final;
Ssubt = Sdemod - Sfit;

% estimate deltaR1
F = zeros(Necho*2, Ncomp);
% Chat(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE);
% Dhat(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE);
% Chat(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE);
% Dhat(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE);
% % !!! Need to make into matrix calculation
% for ee = 1:Necho
%     for nn = 1:Ncomp
%         E1 = exp(-TR*abs(R1_iter(nn)));
%         F(ee,nn) = (rhoR_final(nn)*Chat(ee,nn) - rhoI_final(nn)*Dhat(ee,nn))*sin(FA(ee))*E1*TR/(1-cos(FA(ee))*E1);
%         F(Necho+ee,nn) = (rhoR_final(nn)*Dhat(ee,nn) + rhoI_final(nn)*Chat(ee,nn))*sin(FA(ee))*E1*TR/(1-cos(FA(ee))*E1);
%     end
% end
for nn = 1:Ncomp
    E1 = exp(-TR*abs(R1_iter(nn)));
    F(:,nn) = TR*E1/(1-E1)*([C(:,nn), -D(:,nn); D(:,nn), C(:,nn)]*[rhoR_final(nn);rhoI_final(nn)]);
end

r1 = (F'*F)\F'*Ssubt;

R1_delta = r1;
R1_iter = abs(R1_iter + R1_delta);
R1_final = R1_iter;

% calculate fitted data and residual
C(:,1) = cos(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(1))));
D(:,1) = sin(2*pi*0*TE).*exp(-abs(R2s_final(1))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(1))))./(1-cos(FA)*exp(-TR*abs(R1_final(1))));
C(:,2) = (cos(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(2))))./(1-cos(FA)*exp(-TR*abs(R1_final(2))));
D(:,2) = (sin(2*pi*TE*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE)...
    .*sin(FA).*(1-exp(-TR*abs(R1_final(2))))./(1-cos(FA)*exp(-TR*abs(R1_final(2))));

A = zeros(Necho*2, Ncomp*2);
A(1:Necho,1:2:end) = C;
A(1:Necho,2:2:end) = -D;
A(Necho+1:2*Necho,1:2:end) = D;
A(Necho+1:2*Necho,2:2:end) = C;

% rho_final = (A'*A)\A'*Sdemod;

Sdemod_final = A*rho_final;
S_fit = Sdemod_final(1:Necho) + 1i*Sdemod_final(Necho+1:2*Necho);
S_fit = S_fit.*exp(1i*2*pi*phi_final*TE);

end