function [im, header] = precon_3dute_pfile_bartv300_allec(pfile, ...
    coils, undersamp, ...
    skip, freq_shift, echoes,reg_coe, skip_calib_coil, cc_coil, rNecho,ind_echo_recon, espirit_recon);
% [im, header, rhuser, data, data_grid] = recon_3dute_pfile(pfile,
% coils, undersamp, skip, freq_shift, echoes)
%
% Reconstructs 3D UTE PR image acquired with half-projections and
%   ramp sampling from pfile only.  Reads in scan parameters from pfile.
% INPUTS:
%   pfile - points to a scanner raw data file
%       either use Pfile 'PXXXXX.7', or, for multiple files, can also be
%       'filename01', where the number increases (can use rename_pfiles.x to appropriately rename)
%       use a cell array of the first Pfiles for averaging
%   coils (optional) - can select which coils to reconstruct
%   undersamp (optional) - [undersamp_dc, undersamp_imsize].  Undersampling ratios for
%      density compensation and/or  reconstructed image size
%      (default is [1.0 1.0])
%   skip (optional) - shifts readout by this number of points.  Can
%      be positive or negative
%   freq_shift (optional) - demodulates by this frequency (Hz) to correct for bulk frequency miscalibrations
%   echoes (optional) - can select which echoes to reconstruct
%   rNecho - real number of echo
%   ind_echo_recon  - recon each individual echo sequentially
%   espirit_recon - do espirit recon or just nufft
% OUTPUTS:
%   im - 3d image
%   header - header data from pfile
%   data - raw projection data
%   data_grid - gridded data
%
% Peder Larson 7/28/2008, 6/24/2011
% pcao has some name conflictions for rawloadX_cp (in the src dir) vs. rawloadX and
% rawheadX_cp (in the src dir) vs. rawheadX

if ~isdeployed
    addpath /home/pcao/src/bart-0.2.07/matlab/
    addpath /home/pcao/src/bart-0.3.01/matlab/
    addpath /home/plarson/matlab/3DUTE-recon
%    addpath /netopt/share/ese/ESE_DV26.0_R01/tools/matlab/read_MR/
end

disp(pfile)
header = read_MR_headers(pfile); % <<< DCM >>> 20180312

if (nargin < 2)
    coils = [];
end

if (nargin < 3) || (isempty(undersamp))
    undersamp = 1/header.rdb_hdr.user26;
    undersamp_dc = undersamp; undersamp_imsize = undersamp;
elseif length(undersamp) == 1
    undersamp_dc = undersamp; undersamp_imsize = undersamp;
else
    undersamp_dc = undersamp(1); undersamp_imsize = undersamp(2);
end    

if (nargin < 4) || (isempty(skip))
    skip = 0;
end

if (nargin < 5) || (isempty(freq_shift))
    freq_shift = 0;
end


if (nargin < 6) || (isempty(echoes))
    Necho = header.rdb_hdr.nechoes; % <<< DCM >>> 20180312
    echoes = 1:Necho;
else
    Necho = length(echoes(:));
end


% if isempty(coils)
%   Ncoils = (header.rdb_hdr.dab(2)-header.rdb_hdr.dab(1))+1; % <<< DCM >>> 20180312
%   coils = 1:Ncoils;
% else
     Ncoils = length(coils(:));
% end


if (nargin < 7) || (isempty(reg_coe))
    reg_coe = '-r0.05';
end

if (nargin < 12) || (isempty(espirit_recon))

    espirit_recon = 1;
end


zeropad_factor = 1;

% <<< DCM >>> 20180312 -- KLUDGE ... add all echoes and slices

%frsize = header.rdb_hdr.da_xres; % <<< DCM >>> 20180312
%nframes = header.rdb_hdr.da_yres -1 ; % <<< DCM >>> 20180312 -- possibly off by one because of "baseline" scan
%nslices = header.rdb_hdr.nslices; % <<< DCM >>> 20180312

frsize = header.rdb_hdr.frame_size;%header2.frsize;
nframes = header.rdb_hdr.nframes;%header2.nframes;
nslices = header.rdb_hdr.nslices;%header2.nslices;


nfphases = header.image.fphase;

data = read_MR_rawdata(pfile,'db', 1:nfphases, echoes, 1:nslices,coils); % <<< DCM >>> 20180312 -- replaced rawloadX_cp
data = squeeze(data);
disp('READ DATA SIZE: ')
disp(size(data))

Nramp = header.rdb_hdr.user11;%rhuser(12);
spres = header.rdb_hdr.user1;%rhuser(2);
resz_scale = header.rdb_hdr.user2; %rhuser(3);
FOV = [header.rdb_hdr.user16, header.rdb_hdr.user17, header.rdb_hdr.user18]; %rhuser(17:19).';
Nprojections = header.rdb_hdr.user9;%rhuser(10);
acqs = header.rdb_hdr.user19;%rhuser(20);
shift = [header.rdb_hdr.user20, header.rdb_hdr.user21,0] / spres; %[rhuser(21:22);0].'/spres;
imsize = FOV*10/spres * zeropad_factor / undersamp_imsize;
final_imsize=round(imsize);
a = 1.375;  W = 5; S = calc_kerneldensity(1e-4, a);
gridsize = round(a*imsize);


% For extracting iamges within FOV
rnum = gridsize(1);  cnum = gridsize(2); snum = gridsize(3);
ru_skip = ceil((rnum-final_imsize(1))/2);  rd_skip = rnum - final_imsize(1)-ru_skip;
cu_skip = ceil((cnum-final_imsize(2))/2); cd_skip = cnum - final_imsize(2)-cu_skip;
su_skip = ceil((snum-final_imsize(3))/2); sd_skip = snum - final_imsize(3)-su_skip;

% frsize = size(data,1);
% nframes = size(data,2); 
% Necho = size(data,3);
% nslices = size(data,4);
% Ncoils = size(data,5);

if nframes*nslices >= Nprojections
    % data is now frsize, 2*ceil(Nprojections/(Nslices*2)) echoes,
    % Nslices, Ncoils
    if Nprojections*2 < nframes*nslices
        % storing extra echos in slices
        data = cat(3, data(:,:,:,1:nslices/2,:), data(:,:,:,nslices/2+[1:nslices/2],:));
        nslices = nslices/2;
        Necho = 2*Necho; % num_utes?
        echoes = 1:Necho;
    end
    
    % <<< DCM >>> 20180312 -- reshape and permute the data (interleaved readouts corrected)
    data = permute(data, [2 1 4 3 5]);  %pcao changed coil and echo dimensions
    reordered_projections = [[1:2:nslices],[2:2:nslices]];
    fprintf('DEBUG STUFF: \n');
    disp(size(data));
    %disp(reordered_projections);
    fprintf('frsize = %d, nframes = %d, nslices= %d, Necho = %d, Ncoils = %d; product = %d\n', frsize, nframes, nslices, Necho, Ncoils, frsize*nframes*nslices*Necho*Ncoils);
    fprintf('END DEBUG STUFF\n');
    data = data(:, :, reordered_projections, :);
    data = reshape(data, [frsize nframes*nslices Necho Ncoils]); 
    data = data(:,1:Nprojections,:,:);
        
    if Nprojections*2 < nframes*nslices
        %determine and set the real Number of echo, pcao20170214
        if (nargin < 10) || (isempty(rNecho) || (rNecho < 1))
            sum_nonzero_echo = squeeze(sum(sum(sum(abs(data),1),2),4))> 0;
            Necho = sum(sum_nonzero_echo);
        else
            Necho = rNecho;
        end
        echoes = 1:Necho;
        data = data(:,:,1:Necho,:);
    end
else
    % legacy recon code
    data = read_multiple_pfiles(pfile);
end

% apply different frequency demodulation
if freq_shift ~= 0
    dt = 1/(2*header.rdb_hdr.bw*1e3);
    t = [0:frsize-1].' * dt;
    Sdata = size(data);
    data = data .* repmat( exp(-i*2*pi*freq_shift*t), [1 Sdata(2:end)]);
end

% Determine trajectory
[theta, phi, kmax, dcf] = calc_3dpr_ellipse(FOV*10, spres, spres*resz_scale);

x = cos(phi) .* sin(theta) .* kmax;
y = sin(phi) .* sin(theta) .* kmax;
z = cos(theta) .* kmax;

kscale = 0.5 / max(abs(kmax(:)));
x = kscale * x;
y = kscale * y;
z = kscale * z;

% skip samples?
if skip > 0
    frsize = frsize - skip;
    data = data(1+skip:end,:,:,:);
    % elseif skip < 0
    %   frsize = frsize - skip;
    %   data = [repmat(data(1,:,:,:,:), [-skip 1 1 1 1]); data];
end


[ksp, dcf_all] = calc_pr_ksp_dcf([x(:),y(:),z(:)],Nramp,frsize,dcf,undersamp_dc);
clear x; clear y; clear z;
clear theta; clear phi; clear dcf;
clear kmax;

% ksp((frsize * Nprojections + 1):end,:) = [];
% dcf_all(:,(Nprojections + 1):end) = [];

if skip < 0
    data  = data(1:end+skip,:,:,:);
    dcf_new = dcf_all(1-skip:end,:);
    dcf_new(1,:) = sum(dcf_all(1:1-skip,:),1);  % add density of central skipped points
    dcf_all = dcf_new;
    clear dcf_new
    
    ksp_new = zeros((frsize+skip)*Nprojections,3);
    for n = 1:Nprojections
        ksp_new([1:frsize+skip] + (n-1)*(frsize+skip),:) = ...
            ksp([1-skip:frsize] + (n-1)*frsize,:);
    end
    ksp = ksp_new;
    clear ksp_new
    
    %    frsize = frsize + skip; % not necessary to change
end

if (nargin < 8) || (isempty(skip_calib_coil))
    skip_calib_coil = 0;%skip the coil calibration
end

if (nargin < 9) || (isempty(cc_coil))
    cc_coil = 0; % the coil compression 
end



if cc_coil && (length(coils) > cc_coil) %do coil compression
    disp('Coil compression')
    data(:,:,:,1:cc_coil) = bart(sprintf('cc -r12 -P%d -S',cc_coil), data);
%     clear data;
%     data = cc_data;
%     clear cc_data;
    coils = 1:cc_coil;
    skip_calib_coil = 0;%cannot skip the sensitivity measurement
    data(:,:,:,(cc_coil +1):end) = [];
end


ktraj_all=reshape(ksp,[frsize Nprojections 3]);

ktraj(1,:,:)=ktraj_all(:,:,1)*imsize(1);
ktraj(2,:,:)=ktraj_all(:,:,2)*imsize(2);
ktraj(3,:,:)=ktraj_all(:,:,3)*imsize(3);
tot_npts=frsize*Nprojections;

% dcf_all2(1,:,:)=dcf_all;
% dcf_all2(2,:,:)=dcf_all;
% dcf_all2(3,:,:)=dcf_all;
% clear dcf_all;
% dcf_all = dcf_all2;


for e = 1:length(echoes)
    disp(['Preparing reconstructing echo ' int2str(e) '...'])
    for Ic = 1:length(coils)
        %         disp(['  Reconstructing coil ' int2str(coils(Ic)) '...'])
        %         tic
        
        data_c = squeeze(data(:,:,e,Ic));
        data_pc(:,Ic) = data_c(:).*exp(j*2*pi*(ksp(:,1)*shift(1) + ksp(:,2)*shift(2) + ksp(:,3)*shift(3)));
    end
    clear data_c;
    if e == echoes
        clear data; clear ksp;
    end

    data_pc = conj(data_pc);
    
    tic

    if ~espirit_recon
        disp(['  Reconstructing echo ' int2str(e) '...'])
        im(:,:,:,:,e)=squeeze(bart('nufft -a -p', reshape(dcf_all,[1 tot_npts]), reshape(ktraj, [3 tot_npts]), reshape(data_pc,[1 tot_npts 1 length(coils)])));
    else
        root_dir = pwd;
        list = exist([root_dir '/smap_m1.mat'], 'file');
        if e == 1
            if ~skip_calib_coil || ~list
                disp('nufft to generate calibration k-space')
                %         name = ['/data/larson/brain_uT2/2016-04-27_7T-vounteer/tmpfftec' int2str(echoes)];
                im_under=bart('nufft -a -p', reshape(dcf_all,[1 tot_npts]), reshape(ktraj, [3 tot_npts]), reshape(data_pc,[1 tot_npts 1 length(coils)]));
                k_calb=bart('fft -u 7',im_under);
                k_calb = bart(sprintf('crop 0 %d', 2*round(size(k_calb,1)*0.2)),k_calb);
                k_calb = bart(sprintf('crop 1 %d', 2*round(size(k_calb,2)*0.2)),k_calb);
                k_calb = bart(sprintf('crop 2 %d', 2*round(size(k_calb,3)*0.2)),k_calb);
                k_calb_zerop = padarray(k_calb, round([size(im_under)/2-size(k_calb)/2]));
                clear im_under; clear k_calb;
                
                smap_m1=bart('ecalib -k4 -r12 -m1 -c0.80', k_calb_zerop);  %two sets sensitivity maps are needed here, as tested on the /data/vig2/UTE_ZTE/3dute/brain/20150506_TEphase
                %         figure, imshow3(abs(squeeze(smap_m1(:,:,2,:))));
                clear k_calb_zerop;
                %         if ~isdeployed
                save('smap_m1.mat','smap_m1');
                %         end
            else
                load smap_m1.mat
            end
        end
        
        smapall(:,:,:,:,1,e) = smap_m1;
        dataall(:,:,:,:,1,e) = reshape(data_pc,[1 tot_npts 1 length(coils)]);
        ktrjall(:,:,1,1,1,e) = reshape(ktraj, [3 tot_npts]);
        decfall(:,:,1,1,1,e) = reshape(dcf_all,[1 tot_npts]);
    end
    toc
end

if espirit_recon
    clear data
    disp('Reconstructing all echos ')
    %try add ' -o' scale *; add -n turn off random shift; add -I to choose
    %iterative thresholding; move -l1 to reg_coe
    %     recon_l1 =  bartv207(['nusense -I -o -n ' reg_coe ' -p'], reshape(dcf_all,[1 tot_npts]), reshape(ktraj, [3 tot_npts]), reshape(data_pc,[1 tot_npts 1 length(coils)]),smap_m1);
    bartcmd = ['pics ' reg_coe];
    bartaddcmd{1} = '-p ';
    bartaddcmd{2} = '-t ';
    bartaddcmd{3} = ' ';
    bartaddcmd{4} = ' ';
    
    if nargin < 11 || (isempty(cc_coil))
        ind_echo_recon = 0;
    end
    
    if ind_echo_recon
        for e = 1:length(echoes)
            disp(['  Individual reconstructing echo ' int2str(e) '...'])
            recon_l1(:,:,:,1,1,e) =  bartv301addcmd(bartcmd, bartaddcmd, decfall(:,:,1,1,1,e), ktrjall(:,:,1,1,1,e), dataall(:,:,:,:,1,e),smapall(:,:,:,:,1,e));
        end
    else
        recon_l1 =  bartv301addcmd(bartcmd, bartaddcmd, decfall, ktrjall, dataall,smapall);
    end
    im = squeeze(recon_l1); %somehow recon has two sets of data, like corespond to two sets of sensitivity maps
    
    clear recon_l1
end
return

function data = read_multiple_pfiles(pfile)
% legacy code for reading multiple pfile data

MAX_FRAMES = 16384;

if ~iscell(pfile)
    temp = pfile;
    pfile = cell(1);
    pfile{1} = temp;
end

nex = length(pfile);

[data1, header, rhuser] = rawloadX(pfile{1}, [0:MAX_FRAMES],1,1);
Nprojections = rhuser(10);
acqs = rhuser(20);

frsize = size(data1,1);
Ncoils = size(data1,5);

data = zeros(frsize, Nprojections, Necho, Ncoils);

for n = 1:nex
    for a = 1:acqs
        pfile_a = parse_pfile(pfile{1}, a);
        
        disp(['Reading ' pfile_a '...'])
        tic
        [data1] = rawloadX(pfile_a, [0:MAX_FRAMES],1,1);
        toc
        
        data(:, [MAX_FRAMES*(a-1)+1:min(Nprojections, MAX_FRAMES*a)], :,:) = ...
            data(:, [MAX_FRAMES*(a-1)+1:min(Nprojections, MAX_FRAMES*a)], :,:) + ...
            squeeze(data1(:, 1:min(Nprojections - (a-1)*MAX_FRAMES,MAX_FRAMES),:,1,:));
        clear data1;
        
    end
end

return


function pfile_name = parse_pfile(pfilestring, acq)

% determine if PXXXXX.7 filename or other
if strcmp(pfilestring(end-1:end), '.7')
    pfile_num = sscanf(pfilestring(end-6:end-2),'%d');
    pfile_path = pfilestring(1:end-8);
    pfile_name = sprintf('%sP%05d.7', pfile_path, pfile_num + acq-1);
else
    pfile_num = sscanf(pfilestring(end-1:end),'%d');
    if isempty(pfile_num) % just single pfile
        pfile_name = pfilestring;
    else
        pfile_path = pfilestring(1:end-2);
        pfile_name = sprintf('%s%02d', pfile_path, pfile_num + acq-1);
    end
end

return
