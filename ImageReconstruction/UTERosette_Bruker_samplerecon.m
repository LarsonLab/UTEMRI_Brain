
if (~exist('Gmri','file'))
    %set up IRT 
    currentDirectory = pwd;
    cd('/home/xshen/Documents/MATLAB/irt/irt/'); setup
    cd(currentDirectory)
end

data_directory='/working/larson3/Xin/Preclinical_UTE/8/8/'; %directory with the data to reconstruct

% make sure 'method' file is modified if necessary
method_data=readparam_Bruker(fullfile(data_directory,'method'));
if(contains(method_data.TITLE,'360'))
    fidname='rawdata.job0';
    pv360=1;
else
    fidname='fid';
    pv360=0;
end

acqp=readparam_Bruker(fullfile(data_directory,'acqp'));
fid=getfid_Bruker(fullfile(data_directory,fidname));

if pv360
    nch=method_data.PVM_EncNReceivers;
    ns=method_data.ShotPoints+method_data.Deadpoints+method_data.DeadpointsEnd; 
else
    nch=sum(method_data.PVM_EncActReceivers=='n');
    ns =(acqp.ACQ_size(1)/2)/method_data.Shots;
end

nslices=acqp.NSLICES;
proj=method_data.NPro; 
shots=method_data.Shots; 
coils=nch;

specs=1; 
if(method_data.MEGAOn)
   specs=method_data.MegaListSize;
   if contains(method_data.MegaWSOnly,'Yes')
       %there are two sequence versions: older version forced 2x the number
       %of acquisitions with MEGA of any kind on. New version will just do
       %1 phase encode for each shot if the mega is ws only.

       if (numel(fid)/(proj*nslices))>=(ns*nch*shots*method_data.MegaListSize)
           specs=method_data.MegaListSize;
       else
           specs=1;
       end
   end
end

datablock=numel(fid)/(proj * max([nslices specs]));
fid=reshape(fid, [datablock max([nslices specs]) proj]);
fid=fid(1:(ns*nch*shots),:,:);
fid=reshape(fid, [ns*shots nch max([nslices specs]) proj]);

k=petalutegradsr(method_data,acqp);

sikmax=1./(2*method_data.PVM_SpatResol);

dims=3; slices=1; 
if(contains(method_data.Acq2D,'Yes'))
    k(3,:,:)=0;
    dims=2; 
    slices=method_data.PVM_SPackArrNSlices; 
end

% sets reconstruction size based on type of image 
if shots>1
    Np=64;
elseif dims==2
    Np=512; 
else 
    Np=256;
end
split=round(ns/2)-1;

% controls delay parameters
ein = [1 split-1; split-5 ns-5];
eout= [2 split;   1 1-split+ns];
fidp= [2 split+1;   split ns];
szk = [split ns-split+1];

% ein = [1 split; split-5 ns-5];
% eout= [1 split;   1 1-split+ns];
% fidp= [3 split+2;   split ns];
% szk = [split ns-split+1];

if shots==1 
    if dims==2 
        combined=zeros([Np Np slices 2]); %2 echoes
        mask2 = true([Np Np]);
        nufft2 = {[Np Np], [3 3], 2*[Np Np], [Np/2 Np/2], 'table', 2^10, 'minmax:kb'};
        f.basis = {'rect'};

        recon=zeros([size(mask2) slices coils]);
        
        for rec=1:2 % first and second echoes        

            sz=[szk(rec) size(k,3)];
            kx2=zeros(sz);ky2=zeros(sz);
        
            kx2(eout(rec,1):eout(rec,2),:)=k(1,ein(rec,1):ein(rec,2),:);
            ky2(eout(rec,1):eout(rec,2),:)=k(2,ein(rec,1):ein(rec,2),:);
        
            Gm1 = Gmri([kx2(:) ky2(:)], mask2, 'fov', method_data.PVM_Fov', 'basis', f.basis, 'nufft', nufft2);
            kdens=(ir_mri_density_comp([kx2(:) ky2(:)],'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',method_data.PVM_Fov'}))';
            
            for slice=1:slices %slices
                for coil=1:coils %coil
                    myfid=squeeze(fid(fidp(rec,1):fidp(rec,2),coil,slice,:));
                    xcp  =   (Gm1'*(((myfid(:)).*(kdens)')));
                    recon(:,:,slice,coil) =embed(xcp,mask2);
                end
            end
            
            if coils>1
            for slice=1:slices
                px = angle(squeeze((recon(:,:,slice,:,1))));
                nws_water_nuf=squeeze((recon(:,:,slice,:,1))).* exp( -1i * px );
                ll=estimate_csm_walsh(nws_water_nuf);
                ll(ll < eps) = 1;
                combined(:,:,slice,rec)= sum(conj(ll).*nws_water_nuf ,3);
            end
            else
                combined(:,:,slice)=recon;
            end
        end
        [~,b]=sort(method_data.PVM_ObjOrderList);
        combined(:,:,:,:)=combined(:,:,b,:);
    else % 3D image 
        combined=zeros([Np Np Np 2]); %2 echoes
        mask2 = true([Np Np Np]);
        nufft2 = {[Np Np Np], [3 3 3], 2*[Np Np Np], [Np/2 Np/2 Np/2], 'table', 2^10, 'minmax:kb'};
        f.basis = {'rect'};
        recon=zeros([size(mask2) coils]);
        
        for rec=1:2 % first and second echoes
            sz=[szk(rec) size(k,3)];
            % kx3=zeros(sz);ky3=zeros(sz);kz3=zeros(sz);
            % kx4=zeros(sz);ky4=zeros(sz);kz4=zeros(sz);
            % 
            % kx3(eout(rec,1):eout(rec,2),:)=k(1,ein(rec,1):ein(rec,2),:);
            % ky3(eout(rec,1):eout(rec,2),:)=k(2,ein(rec,1):ein(rec,2),:);
            % kz3(eout(rec,1):eout(rec,2),:)=k(3,ein(rec,1):ein(rec,2),:);
            % kx4(eout(rec,1)-1:eout(rec,2),:)=k(1,ein(rec,1):ein(rec,2)+1,:);
            % ky4(eout(rec,1)-1:eout(rec,2),:)=k(2,ein(rec,1):ein(rec,2)+1,:);
            % kz4(eout(rec,1)-1:eout(rec,2),:)=k(3,ein(rec,1):ein(rec,2)+1,:);
            % kx2=(kx3+kx4)/2;
            % ky2=(ky3+ky4)/2;
            % kz2=(kz3+kz4)/2;
            % kx2=pi*kx2/max(kx2(:));
            % ky2=pi*ky2/max(ky2(:));
            % kz2=pi*kz2/max(kz2(:));
            kx2=zeros(sz);ky2=zeros(sz);kz2=zeros(sz);
        
            kx2(eout(rec,1):eout(rec,2),:)=k(1,ein(rec,1):ein(rec,2),:);
            ky2(eout(rec,1):eout(rec,2),:)=k(2,ein(rec,1):ein(rec,2),:);
            kz2(eout(rec,1):eout(rec,2),:)=k(3,ein(rec,1):ein(rec,2),:);
        
            Gm1 = Gmri([kx2(:) ky2(:) kz2(:)], mask2, 'fov', method_data.PVM_Fov', 'basis', f.basis, 'nufft', nufft2);
            kdens=(ir_mri_density_comp([kx2(:) ky2(:) kz2(:)],'pipe','G',Gm1.arg.Gnufft,'arg_pipe',{'fov',method_data.PVM_Fov'}))';
            
            for coil=1:coils %coil
                myfid=squeeze(fid(fidp(rec,1):fidp(rec,2),coil,1,:,:));
                xcp  =   (Gm1'*(((myfid(:)).*(kdens)')));
                recon(:,:,:,coil) =embed(xcp,mask2);
            end

            if coils>1
                px = angle(squeeze((recon(:,:,:,:,1))));
                nws_water_nuf=squeeze((recon(:,:,:,:,1))).* exp( -1i * px );
                ll=estimate_csm_walsh(nws_water_nuf);
                ll(ll < eps) = 1;
                combined(:,:,:,rec)= sum(conj(ll).*squeeze(nws_water_nuf),4);
            else
                combined(:,:,:,rec)= recon;
            end
        end
    end
end
save('image_8.mat','combined')