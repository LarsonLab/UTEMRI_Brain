function [k,t,ktraj,g,gtraj,theoryktraj,theoryk,kgradtime,gradtime,kinttraj]=petalutegrads(meth,acqp)
%t is rftime 

t=[];k=[];ktraj=[];g=[];gtraj=[];
if(contains(meth.TITLE,'360'))
    pv360=1;
    rfpoints=meth.ShotPoints+meth.Deadpoints+meth.DeadpointsEnd;
else
    rfpoints=(acqp.ACQ_size(1)/2)/meth.Shots; 
    fidname='fid';
end

npetals=meth.NPro; %acqp.ACQ_size(2);


try
    shots=meth.Shots;
catch
    shots=1; 
end

try
    deend=meth.DeadtimeEnd;
    dpend=meth.DeadpointsEnd;
catch
    deend=0;
    dpend=0; 
end

acqdim=3;
try 
    if contains(meth.Acq2D,'Yes')
        acqdim=2; 
    end
catch
end

%gradres=(meth.PVM_AcquisitionTime+meth.Deadtime)/numel(meth.GradShapec1);
%+meth.Deadtime+deend remove this from time in the below line 
gradtime=linspace(0,meth.PVM_AcquisitionTime, numel(meth.GradShapec1));%(1:(numel(meth.GradShapec1)))*gradres; %ms
%rftime=  linspace(0,meth.PVM_AcquisitionTime+meth.Deadtime+deend, meth.Deadpoints+meth.PVM_Matrix(1)+dpend); 
rftime= 1000*(1:(rfpoints))/meth.PVM_EffSWh;    %ms - modified 5/9/2023 call with uzay 
t=rftime; 

kmax=meth.PVM_Matrix./(2*meth.PVM_Fov);

%interpolate gradient time to RF time  
rfgs1c2=interp1(gradtime,meth.GradShapes1c2,rftime,'linear',0);
rfgs1s2=interp1(gradtime,meth.GradShapes1s2,rftime,'linear',0);
rfgc1c2=interp1(gradtime,meth.GradShapec1c2,rftime,'linear',0);
rfgc1s2=interp1(gradtime,meth.GradShapec1s2,rftime,'linear',0);
rfgc1=  interp1(gradtime,meth.GradShapec1,  rftime,'linear',0);

%integrate to get k routes : factors of 1e3 for ms
rfks1c2=cumtrapz(rftime/1e3, rfgs1c2);
rfks1s2=cumtrapz(rftime/1e3, rfgs1s2);
rfkc1c2=cumtrapz(rftime/1e3, rfgc1c2);
rfkc1s2=cumtrapz(rftime/1e3, rfgc1s2);
rfkc1=  cumtrapz(rftime/1e3, rfgc1);

%combine for each petal 
% if nargout>1
    k=zeros([acqdim numel(rftime) npetals]); % 1/mm
    
    for p=1:npetals
        k(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(meth.GradAmpR1(p)*rfkc1c2-meth.GradAmpR2(p)*rfkc1s2-meth.GradAmpP1(p)*rfks1s2-meth.GradAmpP2(p)*rfks1c2);
        k(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(meth.GradAmpR2(p)*rfkc1c2+meth.GradAmpR1(p)*rfkc1s2-meth.GradAmpP2(p)*rfks1s2+meth.GradAmpP1(p)*rfks1c2);
        if acqdim>2
            k(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(meth.GradAmpS(p)*rfkc1);
        end
    end

    k=repmat(k,[1 shots 1]);
% end

% if nargout>2
% 
%     ntraj=meth.TrajBeta*meth.TrajPhi;
%     trajR1=zeros([1 ntraj]);trajR2=zeros([1 ntraj]);trajP1=zeros([1 ntraj]);trajP2=zeros([1 ntraj]);trajS=zeros([1 ntraj]);
% 
%     n=1;
%     nrm=sqrt(meth.Petalw1.^2+meth.Petalw2.^2);
% 
%     if ~isfield(meth,'TrajPhiVals') %this was a later addition to the method
%         for i=1:meth.TrajPhi
%             cosp=cos(-pi/2 + (i-1)*pi/meth.TrajPhi);
%             sinp=sin(-pi/2 + (i-1)*pi/meth.TrajPhi);
%             for j=1:meth.TrajBeta
%                 cosb=cos((j-1)*2*pi/meth.TrajBeta);
%                 sinb=sin((j-1)*2*pi/meth.TrajBeta);
%                 trajR1(n)=cosp*cosb*meth.Petalw1/nrm; 
%                 trajR2(n)=cosp*sinb*meth.Petalw1/nrm; 
%                 trajP1(n)=cosp*cosb*meth.Petalw2/nrm; 
%                 trajP2(n)=cosp*sinb*meth.Petalw2/nrm; 
%                 trajS(n) =sinp; 
%                 n=n+1;
%             end
%         end
%     else
%         for i=1:meth.TrajPhi
%             cosp=cos(meth.TrajPhiVals(i));
%             sinp=sin(meth.TrajPhiVals(i));
%             for j=1:meth.TrajBeta
%                 cosb=cos(meth.TrajBetaVals(j));
%                 sinb=sin(meth.TrajBetaVals(j));
%                 trajR1(n)=cosp*cosb*meth.Petalw1/nrm; 
%                 trajR2(n)=cosp*sinb*meth.Petalw1/nrm; 
%                 trajP1(n)=cosp*cosb*meth.Petalw2/nrm; 
%                 trajP2(n)=cosp*sinb*meth.Petalw2/nrm; 
%                 trajS(n) =sinp; 
%                 n=n+1;
%             end
%         end
%     end
% 
%     ktraj=zeros([acqdim numel(rftime) ntraj]); % 1/mm
% 
%     for p=1:ntraj
%         ktraj(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(trajR1(p)*rfkc1c2-trajR2(p)*rfkc1s2-trajP1(p)*rfks1s2-trajP2(p)*rfks1c2);
%         ktraj(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(trajR2(p)*rfkc1c2+trajR1(p)*rfkc1s2-trajP2(p)*rfks1s2+trajP1(p)*rfks1c2);
%         if acqdim>2
%             ktraj(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(trajS(p) *rfkc1);
%         end
%     end
% end
% 
% if nargout>3
%     g=zeros([acqdim numel(rftime) npetals]); %Hz/mm
% 
%     for p=1:npetals
%         g(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(meth.GradAmpR1(p)*rfgc1c2-meth.GradAmpR2(p)*rfgc1s2-meth.GradAmpP1(p)*rfgs1s2-meth.GradAmpP2(p)*rfgs1c2);
%         g(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(meth.GradAmpR2(p)*rfgc1c2+meth.GradAmpR1(p)*rfgc1s2-meth.GradAmpP2(p)*rfgs1s2+meth.GradAmpP1(p)*rfgs1c2);
%         if acqdim>2
%             g(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(meth.GradAmpS(p)*rfgc1);
%         end
%     end
% 
%     g=repmat(g,[1 shots 1]);
% end
% 
% 
% if nargout>4
%     gtraj=zeros([acqdim numel(rftime) ntraj]); %Hz/mm
% 
%     for p=1:ntraj
%         gtraj(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(trajR1(p)*rfgc1c2-trajR2(p)*rfgc1s2-trajP1(p)*rfgs1s2-trajP2(p)*rfgs1c2);
%         gtraj(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(trajR2(p)*rfgc1c2+trajR1(p)*rfgc1s2-trajP2(p)*rfgs1s2+trajP1(p)*rfgs1c2);
%         if acqdim>2
%             gtraj(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(trajS(p)*rfgc1);
%         end
%     end
% end
% 
% if nargout>5
% 
%     rf1shot=rftime;%(1:meth.PVM_EncMatrix(1));
% 
%     theoryktraj=zeros([acqdim numel(rf1shot) ntraj]);
%     if ~isfield(meth,'TrajPhiVals')
%         [beta,phi]=ndgrid(((0:meth.TrajBeta-1)*2*pi)/meth.TrajBeta, -pi/2+((0:meth.TrajPhi-1)*pi)/meth.TrajPhi);
%     else
%         beta=meth.TrajBetaVals;
%         phi=meth.TrajPhiVals;
%         [beta,phi]=ndgrid(beta,phi);
%     end
% 
%     for i=1:ntraj
%         kxy=cos(phi(i)).*sin(meth.Petalw1*rf1shot/1000).*exp(1i*(meth.Petalw2*rf1shot/1000+beta(i)));
%         theoryktraj(1,:,i)=kmax(1)*real(kxy);
%         theoryktraj(2,:,i)=kmax(2)*imag(kxy);
%         if acqdim>2
%             theoryktraj(3,:,i)=kmax(3)*sin(phi(i)).*sin(meth.Petalw1*rf1shot/1000);
%         end
%     end
% 
%     theoryktraj=repmat(theoryktraj,[1 shots 1]);
% 
% end
% 
% if nargout>6
% 
%     rf1shot=rftime;%(1:meth.PVM_EncMatrix(1));
% 
%     theoryk=zeros([acqdim numel(rf1shot) npetals]);
%     [beta,phi]=ndgrid(((0:meth.NBeta-1)*2*pi)/meth.NBeta, -pi/2+((0:meth.NPhi-1)*pi)/meth.NPhi);
% 
%     for i=1:npetals
%         kxy=cos(phi(i)).*sin(meth.Petalw1*rf1shot/1000).*exp(1i*(meth.Petalw2*rf1shot/1000+beta(i)));
%         theoryk(1,:,i)=kmax(1)*real(kxy);
%         theoryk(2,:,i)=kmax(2)*imag(kxy);
%         if acqdim>2
%             theoryk(3,:,i)=kmax(3)*sin(phi(i))*sin(meth.Petalw1*rf1shot/1000);
%         end
%     end
% 
%     theoryk=repmat(theoryk,[1 shots 1]);
% end
% 
% if nargout>7
% 
%     %integrate to get k routes : factors of 1e3 for ms
%     gs1c2=cumtrapz(gradtime/1e3, repmat(meth.GradShapes1c2,[1]));
%     gs1s2=cumtrapz(gradtime/1e3, repmat(meth.GradShapes1s2,[1]));
%     gc1c2=cumtrapz(gradtime/1e3, repmat(meth.GradShapec1c2,[1]));
%     gc1s2=cumtrapz(gradtime/1e3, repmat(meth.GradShapec1s2,[1]));
%     gc1=  cumtrapz(gradtime/1e3, repmat(meth.GradShapec1,[1]));
% 
%     %combine for each petal 
% 
%     kgradtime=zeros([acqdim numel(gradtime) npetals]); % 1/mm
% 
%     for p=1:npetals
%         kgradtime(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(meth.GradAmpR1(p)*gc1c2-meth.GradAmpR2(p)*gc1s2-meth.GradAmpP1(p)*gs1s2-meth.GradAmpP2(p)*gs1c2);
%         kgradtime(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(meth.GradAmpR2(p)*gc1c2+meth.GradAmpR1(p)*gc1s2-meth.GradAmpP2(p)*gs1s2+meth.GradAmpP1(p)*gs1c2);
%         if acqdim>2
%             kgradtime(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(meth.GradAmpS(p)*gc1);
%         end
%     end
% 
%     kgradtime=repmat(kgradtime,[1 shots 1]);
% 
% end 
% 
% if nargout>9
%     kinttraj=zeros([acqdim numel(rftime) ntraj]); %Hz/mm
% 
%     for p=1:ntraj
%         kinttraj(1,:,p)=meth.ReadGrad* meth.PVM_GradCalConst/100*(trajR1(p)*rfkc1c2-trajR2(p)*rfkc1s2-trajP1(p)*rfks1s2-trajP2(p)*rfks1c2);
%         kinttraj(2,:,p)=meth.PhaseGrad*meth.PVM_GradCalConst/100*(trajR2(p)*rfkc1c2+trajR1(p)*rfkc1s2-trajP2(p)*rfks1s2+trajP1(p)*rfks1c2);
%         if acqdim>2
%             kinttraj(3,:,p)=meth.SliceGrad*meth.PVM_GradCalConst/100*(trajS(p)*rfkc1);
%         end
%     end
% 
% end