function [Sdemod, Ssubt, S_fit, phi_final, R2s_final, rho_final, phi_delta, R2s_delta] = ...
    IDEALiter_2comp_r2s_multipeak(S0, TE, phi_iter, R2s_iter, multipeak, Necho)

% default Necho is 8
if nargin < 6
    Necho = 8;
end

% determine if multi FA acquisition or not
if length(TE) > Necho
    Nflipangle = length(TE)/Necho;
    if mod(length(TE), Necho) ~= 0
        error('Wrong input of Necho.');
    end
elseif length(TE) == Necho
    Nflipangle = 1;
else
    error('Wrong input of Necho.');
end
% fprintf('Necho: %i, Nflipangle: %i \n', Necho, Nflipangle);

Ncomp = 2;

% rearrange S0 and TE into acquisitions
S0 = reshape(S0,[Necho, Nflipangle]); % Necho x Nflipangle Complex
TE = reshape(TE,[Necho, Nflipangle]); % Necho x Nflipangle

% demodulating B0 field
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE); % Necho x Nflipangle Complex
Sdemod = [real(Sdemod_C); imag(Sdemod_C)]; % Necho x Nflipangle*2 Real
Sdemod = reshape(Sdemod,[Necho*Nflipangle*2,1]); % Necho*Nflipangle*2 Real

% estmate rho given B0
C = cell(1,Nflipangle);
D = cell(1,Nflipangle);
A = cell(1,Nflipangle);
Aall = zeros(Necho*2*Nflipangle, Ncomp*2*Nflipangle);
for FAnumber = 1:Nflipangle
    C{FAnumber}(:,1) = cos(2*pi*0*TE(:,FAnumber)).*exp(-abs(R2s_iter(1))*TE(:,FAnumber));
    D{FAnumber}(:,1) = sin(2*pi*0*TE(:,FAnumber)).*exp(-abs(R2s_iter(1))*TE(:,FAnumber));
    C{FAnumber}(:,2) = (cos(2*pi*TE(:,FAnumber)*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE(:,FAnumber));
    D{FAnumber}(:,2) = (sin(2*pi*TE(:,FAnumber)*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_iter(2))*TE(:,FAnumber));
    
    A{FAnumber} = zeros(Necho*2, Ncomp*2);
    A{FAnumber}(1:Necho,1:2:end) = C{FAnumber};
    A{FAnumber}(1:Necho,2:2:end) = -D{FAnumber};
    A{FAnumber}(Necho+1:2*Necho,1:2:end) = D{FAnumber};
    A{FAnumber}(Necho+1:2*Necho,2:2:end) = C{FAnumber};
    
    Aall((FAnumber-1)*Necho*2+1:FAnumber*Necho*2, (FAnumber-1)*Ncomp*2+1:FAnumber*Ncomp*2) = A{FAnumber};
end

rho = (Aall'*Aall)\Aall'*Sdemod;
rho = reshape(rho,[Ncomp*2, Nflipangle]);
rhoR = rho(1:2:end,:);
rhoI = rho(2:2:end,:);

% calculate residual: subtract estimated signal from acquired signal
Sreal = zeros(size(S0));
Simag = zeros(size(S0));
for FAnumber = 1:Nflipangle
    Sreal(:,FAnumber) = diag(repmat([rhoR(:,FAnumber);rhoI(:,FAnumber)]',[Necho,1])*[C{FAnumber}';-D{FAnumber}']);
    Simag(:,FAnumber) = diag(repmat([rhoR(:,FAnumber);rhoI(:,FAnumber)]',[Necho,1])*[D{FAnumber}';C{FAnumber}']);
end
Sfit = [Sreal; Simag];
Sfit = reshape(Sfit,[Necho*Nflipangle*2,1]);
Ssubt = Sdemod - Sfit;

% estimate delta phi and delta rho
GR = zeros(size(S0));
GI = zeros(size(S0));
for FAnumber = 1:Nflipangle
    GR(:,FAnumber) = 2*pi*TE(:,FAnumber).*diag(repmat([rhoR(:,FAnumber);rhoI(:,FAnumber)]',[Necho,1])*[-D{FAnumber}';-C{FAnumber}']);
    GI(:,FAnumber) = 2*pi*TE(:,FAnumber).*diag(repmat([rhoR(:,FAnumber);rhoI(:,FAnumber)]',[Necho,1])*[C{FAnumber}';-D{FAnumber}']);
end
G = [GR; GI];
G = reshape(G,[Necho*Nflipangle*2,1]);

Ball = [G Aall];

y = (Ball'*Ball)\Ball'*Ssubt;

phi_delta = y(1);
rhoR_delta = y(2:2:end);
rhoR_delta = reshape(rhoR_delta,[Ncomp, Nflipangle]);
rhoI_delta = y(3:2:end);
rhoI_delta = reshape(rhoI_delta,[Ncomp, Nflipangle]);

phi_iter = phi_iter + phi_delta;
phi_final = phi_iter;
rhoR_final = rhoR + rhoR_delta;
rhoI_final = rhoI + rhoI_delta;
rho_final = [rhoR; rhoI];
rho_final(1:2:end,:) = rhoR_final;
rho_final(2:2:end,:) = rhoI_final;

% demodulating B0 field again
Sdemod_C = S0.*exp(-1i*2*pi*phi_iter*TE);
Sdemod = [real(Sdemod_C); imag(Sdemod_C)]; % Necho x Nflipangle*2 Real
Sdemod = reshape(Sdemod,[Necho*Nflipangle*2,1]); % Necho*Nflipangle*2 Real

% calculate residual: subtract estimated signal from acquired signal
Sreal = zeros(size(S0));
Simag = zeros(size(S0));
for FAnumber = 1:Nflipangle
    Sreal(:,FAnumber) = diag(repmat([rhoR_final(:,FAnumber);rhoI_final(:,FAnumber)]',[Necho,1])*[C{FAnumber}';-D{FAnumber}']);
    Simag(:,FAnumber) = diag(repmat([rhoR_final(:,FAnumber);rhoI_final(:,FAnumber)]',[Necho,1])*[D{FAnumber}';C{FAnumber}']);
end
Sfit = [Sreal; Simag];
Sfit = reshape(Sfit,[Necho*Nflipangle*2,1]);
Ssubt = Sdemod - Sfit;

% estimate delta R2s
Eall = zeros(Necho*2*Nflipangle, Ncomp);
for FAnumber = 1:Nflipangle
    E = zeros(Necho*2, Ncomp);
    rhoR_temp = rhoR_final(:,FAnumber);
    rhoI_temp = rhoI_final(:,FAnumber);
    TE_temp = TE(:,FAnumber);
    % !!! Need to make into matrix calculation
    for ee = 1:Necho
        for nn = 1:Ncomp
            E(ee,nn) = -rhoR_temp(nn)*C{FAnumber}(ee,nn)*TE_temp(ee) + rhoI_temp(nn)*D{FAnumber}(ee,nn)*TE_temp(ee);
            E(Necho+ee,nn) = -rhoR_temp(nn)*D{FAnumber}(ee,nn)*TE_temp(ee) - rhoI_temp(nn)*C{FAnumber}(ee,nn)*TE_temp(ee);
        end
    end
    Eall((FAnumber-1)*Necho*2+1:FAnumber*Necho*2,:) = E;
end

z = (Eall'*Eall)\Eall'*Ssubt;

R2s_delta = z;
R2s_iter = abs(R2s_iter + R2s_delta);
R2s_final = R2s_iter;

% calculate fitted signal and residual
C = cell(1,Nflipangle);
D = cell(1,Nflipangle);
A = cell(1,Nflipangle);
Aall = zeros(Necho*2*Nflipangle, Ncomp*2*Nflipangle);
for FAnumber = 1:Nflipangle
    C{FAnumber}(:,1) = cos(2*pi*0*TE(:,FAnumber)).*exp(-abs(R2s_final(1))*TE(:,FAnumber));
    D{FAnumber}(:,1) = sin(2*pi*0*TE(:,FAnumber)).*exp(-abs(R2s_final(1))*TE(:,FAnumber));
    C{FAnumber}(:,2) = (cos(2*pi*TE(:,FAnumber)*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE(:,FAnumber));
    D{FAnumber}(:,2) = (sin(2*pi*TE(:,FAnumber)*multipeak.chemshift)*multipeak.alpha').*exp(-abs(R2s_final(2))*TE(:,FAnumber));
    
    A{FAnumber} = zeros(Necho*2, Ncomp*2);
    A{FAnumber}(1:Necho,1:2:end) = C{FAnumber};
    A{FAnumber}(1:Necho,2:2:end) = -D{FAnumber};
    A{FAnumber}(Necho+1:2*Necho,1:2:end) = D{FAnumber};
    A{FAnumber}(Necho+1:2*Necho,2:2:end) = C{FAnumber};
    
    Aall((FAnumber-1)*Necho*2+1:FAnumber*Necho*2, (FAnumber-1)*Ncomp*2+1:FAnumber*Ncomp*2) = A{FAnumber};
end

rho_final = reshape(rho_final,[Ncomp*Nflipangle*2,1]);
Sdemod_final = Aall*rho_final;
Sdemod_final = reshape(Sdemod_final,[Necho*2,Nflipangle]);
S_fit = zeros(size(S0));
for FAnumber = 1:Nflipangle
    Sdemod_temp = Sdemod_final(:,FAnumber);
    S_fit_temp = Sdemod_temp(1:Necho) + 1i*Sdemod_temp(Necho+1:2*Necho);
    S_fit_temp = S_fit_temp.*exp(1i*2*pi*phi_final*TE(:,FAnumber));
    S_fit(:,FAnumber) = S_fit_temp;
end

S_fit = reshape(S_fit,[Necho*Nflipangle,1]);

end