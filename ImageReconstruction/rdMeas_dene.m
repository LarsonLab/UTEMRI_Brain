function [mrprot, mdh, fid] = rdMeas(varargin)
% read in raw fid data from Siemens meas.out/meas.asc file
% returns MrProt, YAPS buffer, MDH entries, and raw fid data
%   E. Auerbach, CMRR, 2016
% usage: [mrprot, mdh, fid] = rdMeas(outPath, ascPath OR multi-RAID idx, reInterpolate=true, deOversample=true, doChannelScaling=false, doFlip=true, doRoFT=false, readTrcRef=false)
% dependencies: parse_mrprot.m
%               parse_xprot.m
%               catstruct.m
%               c_str.m
%
% Example for VA25 (meas.out, meas.asc):
%		[mrprot, mdh, fid] = rdMeas('meas.out', 'meas.asc');
%
% Example for VB15/17/19 (meas.dat):
%		[mrprot, mdh, fid] = rdMeas('meas.dat');
%
% Example for VD11/13 multi-RAID (meas.dat), reading first measurement
%	(default is to read only the last measurement):
%		[mrprot, mdh, fid] = rdMeas('meas.dat', 1);
%

tstart = tic;

version = '2016.11.04';

% ////////////////////////////////////////////////////
% ////// first, some options
% ////////////////////////////////////////////////////

% defaults
outPath = 'meas.out';
ascPath = 'meas.asc';
mraid_idx = -1;

% EPI is often sampled on the ramps and requires interpolation
%   set reInterpolate if this should be done as it is read
reInterpolate = true;

% usually, data is oversampled 2X by default
% set deOversample if rdMeas.m should de-oversample the data
%deOversample = false;
deOversample = false;

% set doChannelScaling = true to scale data by the FFT_SCALE factor for each
% individual channel
doChannelScaling = 0;

% set doFlip = true to reverse the appropriate EPI lines automatically
doFlip = true;

% set doRoFT = true to do RoFT (or not do inv FT after deoversampling...)
doRoFT = false;

% set readTrcRef = true to read in the MB trace reference lines only (normally skipped)
readTrcRef = false;

% 0 = old straight linear interpolation
% 1 = Updated Readout filter from Steen + linear interpolation
% 2 = spline interpolation
interp_scheme = 2;

% read passed parameters
if (nargin >= 1)
    outPath = char(varargin{1});
    tpos = strfind(lower(outPath),'meas.out');
    if (tpos)
        ascPath = outPath;
        ascPath(tpos+5:tpos+7) = 'asc';
    end
    if (nargin >= 2)
    	if (isnumeric(varargin{2})), mraid_idx = varargin{2}; end
        ascPath = char(varargin{2});
        if (nargin >= 3)
            reInterpolate = varargin{3};
            if (nargin >= 4)
                deOversample = varargin{4};
                if (nargin >= 5)
                    doChannelScaling = varargin{5};
                    if (nargin >= 6)
                        doFlip = varargin{6};
                        if (nargin >= 7)
                            doRoFT = varargin{7};
                            if (nargin >= 8)
                                readTrcRef = varargin{8};
                            end
                        end
                    end
                end
            end
        end
    end
end

% end options -------------------------

% first, open the .out/.dat file and check the header

fprintf('\nrdMeas.m version %s by eja\n--------------------------------------------------------------------------------\n', version)
fprintf('Open %s\n',outPath)

fp = fopen(outPath, 'r', 'ieee-le');

% find file size
fseek(fp, 0, 'eof');
fsize = ftell(fp);
fprintf('File size: %.2f MB\n',fsize/(1024*1024));

fseek(fp, 0, 'bof');
hdrsize = fread(fp, 1, 'uint32');

% assume VB13 or newer to start with
VB13 = true;
VD11 = false;

% if the offset is 32 bytes, this is old (pre-VB13) data, so check for the
% meas.asc file--otherwise, the meas.asc equivalent data is embedded in the
% .dat

if (hdrsize == 32)
    % this must be pre-VB13 data
    % for pre-VB13, skip 32-byte header in meas.out and read in parameter data from meas.asc file (MrProt & YAPS data)
    VB13 = false;
    fprintf('This appears to be pre-VB13 meas.out/meas.asc data\n');
    fprintf('Open %s\n',ascPath);
    fp2 = fopen(ascPath, 'r');
    fprintf('Read MrProt & YAPS\n');
    mparr = fread(fp2, inf, 'uint8=>char');
    fprintf('Close meas.asc\n')
    fclose(fp2);
    mrprot = parse_mrprot(mparr);
elseif (hdrsize == 0)
    % this is probably >=VD11 data (multi-RAID)
    fprintf('This appears to be VD11 or newer multi-RAID data\n');
    VD11 = true;
    nMeas = fread(fp, 1, 'uint32');  % number of measurements in this multi-raid file (<=64)
    % read the MrParcRaidFileEntry structures (152 byte length, always 64 present regardless of whether they are used)
    MeasID = zeros(nMeas,1,'uint32');
    FileID = zeros(nMeas,1,'uint32');
    MeasOffset = zeros(nMeas,1,'uint64');
    MeasLen = zeros(nMeas,1,'uint64');
    PatName = cell(nMeas);
    ProtName = cell(nMeas);
    for x=1:nMeas
        MeasID(x) = fread(fp,1,'*uint32');
        FileID(x) = fread(fp,1,'*uint32');
        MeasOffset(x) = fread(fp,1,'*uint64');
        MeasLen(x) = fread(fp,1,'*uint64');
        PatName{x} = fread(fp,64,'char');
        ProtName{x} = fread(fp,64,'char');
        fprintf('Meas %d: MeasID = %d; FileID = %d; Protocol: %s\n', x, MeasID(x), FileID(x), char(ProtName{x}));
    end
    if (nMeas > 1)
    	if (mraid_idx > 0)
    		fprintf('Multi-RAID: reading meas %d (user request)\n',mraid_idx);
    	else
	    	mraid_idx = nMeas;
	    	fprintf('Multi-RAID: no meas index specified -- reading last meas %d by default\n', mraid_idx);
	    end
    else
    	mraid_idx = 1;
    end
    if ((mraid_idx < 1) || (mraid_idx > nMeas))
	    error('Multi-RAID: meas %d out of range',mraid_idx);
    end
    fseek(fp,double(MeasOffset(mraid_idx)),'bof'); % 8+(152*64)=9736, +504 to align to 512byte (MRPARCRAID_SECT_ALIGN)
    hdrsize = double(fread(fp,1,'uint32')) + double(MeasOffset(mraid_idx));
end

if (VB13)
    % this must be >=VB13 data
    nEvp = fread(fp, 1, 'uint32');  % number of embedded evp files (?)
    % read all of the embedded evp files
    EvpName = cell(nEvp);
    EvpDat = cell(nEvp);
    MeasYapsIdx = 0;
    MeasIdx = 0;
    PhoenixIdx = 0;
    for x=1:nEvp
        EvpName{x} = read_cstr(fp);
        if (strcmp(char(EvpName{x}), 'MeasYaps')), MeasYapsIdx = x; end
        if (strcmp(char(EvpName{x}), 'Meas')), MeasIdx = x; end
        if (strcmp(char(EvpName{x}), 'Phoenix')), PhoenixIdx = x; end
        dsize = fread(fp, 1, 'uint32');
        EvpDat{x} = fread(fp, dsize, 'uint8=>char');
    end
    
    if (MeasIdx == 0) % this is VB13 data
        % for VB13, just find the meas.asc part and parse the mrprotocol
        fprintf('This appears to be VB13 meas.dat data\n');
        if (MeasYapsIdx == 0)
            error('meas.asc data not found within meas.dat!')
        end
        %fprintf('%s',char(EvpDat{MeasYapsIdx}))
        mrprot = parse_mrprot(char(EvpDat{MeasYapsIdx}));
    else % this is >=VB15 data
        % for VB15, we have to parse the XProtocol for YAPS parameters, since they are gone from
        % the text Phoenix protocol, which now only contains the bare minimum of parameters
        % needed for Phoenix.
        % we will still parse the text Phoenix protocol, however, since it
        % contains useful arrays that parse_xprot can not handle yet
        if (~VD11), fprintf('This appears to be VB15 or newer meas.dat data\n'); end
        fprintf('** Parsing Phoenix text protocol\n');
        if (PhoenixIdx == 0)
            error('Phoenix data not found within meas.dat!')
        end
        mrprot = parse_mrprot(char(EvpDat{PhoenixIdx}));
        
        fprintf('** Parsing Meas XProtocol\n');
        mrprot.XProtocol = parse_xprot(char(EvpDat{MeasIdx}'), false);
        
        % copy XProtocol MEAS and YAPS parameters to top level for backward compatibility
        mrprot = catstruct(mrprot.XProtocol.MEAS, mrprot.XProtocol.YAPS, mrprot);
    end
end

fclose(fp);

% determine the maximum number of potential channels (for sanity checks)
% and patch if necessary for VBxx UHF 64RX mod
if ((VB13) && (~VD11))
    if (isfield(mrprot.XProtocol.CONFIG.COMMON.CONTROL, 'MODE_64RX'))
        if (mrprot.XProtocol.CONFIG.COMMON.CONTROL.MODE_64RX >= 1)
            fprintf('This is special VBxx 64RX data! Patching receiver channel properties...\n');

            AdditionalRCCSFactor = 2; % 2X RCCS for 64RX
            nChaOrig = mrprot.iMaxNoOfRxChannels;
            nChaMeas = AdditionalRCCSFactor * nChaOrig;
            mrprot.iMaxNoOfRxChannels = nChaMeas;
            mrprot.XProtocol.YAPS.iMaxNoOfRxChannels = nChaMeas;
            mrprot.XProtocol.YAPS.aiMaxNoOfRxChannelsOfCS(1) = AdditionalRCCSFactor * mrprot.XProtocol.YAPS.aiMaxNoOfRxChannelsOfCS(1);
            mrprot.sProtConsistencyInfo.lMaximumNofRxReceiverChannels = AdditionalRCCSFactor * mrprot.sProtConsistencyInfo.lMaximumNofRxReceiverChannels;
            
            % important parameters are:
            %   mrprot.asCoilSelectMeas.asList(:)
            %   mrprot.asCoilSelectMeas.aFFT_SCALE(:)
            % we won't bother with YAPS.alPrimaryPATModes
            
            if (length(mrprot.asCoilSelectMeas) ~= 1), error('ERROR: Detected multiple coil selects; this case is not handled, aborting...'); end

            % copy parameters from first 32ch to second 32ch
            mrprot.asCoilSelectMeas.asList = [mrprot.asCoilSelectMeas.asList ; mrprot.asCoilSelectMeas.asList];
            mrprot.asCoilSelectMeas.aFFT_SCALE = [mrprot.asCoilSelectMeas.aFFT_SCALE ; mrprot.asCoilSelectMeas.aFFT_SCALE];

            % modify
            for idx = nChaOrig+1:nChaMeas
                mrprot.asCoilSelectMeas.asList(idx).lRxChannelConnected = mrprot.asCoilSelectMeas.asList(idx).lRxChannelConnected + (nChaMeas - nChaOrig);
                mrprot.asCoilSelectMeas.aFFT_SCALE(idx).lRxChannel = mrprot.asCoilSelectMeas.aFFT_SCALE(idx).lRxChannel + (nChaMeas - nChaOrig);
                if (~isempty(mrprot.asCoilSelectMeas.asList(idx).sCoilElementID.tElement))
                    mrprot.asCoilSelectMeas.asList(idx).sCoilElementID.tElement = [mrprot.asCoilSelectMeas.asList(idx).sCoilElementID.tElement 'X'];
                end
            end
        end
    end
end
maxchannels = mrprot.iMaxNoOfRxChannels;
    
% call the mex to read in the MDH
tstartmdh = tic;
[mdh, mdh_trc] = rdMeas_mex_mdh(outPath, hdrsize, maxchannels, VD11, VB13, readTrcRef);
telapsedmdh = toc(tstartmdh);
fprintf('   ** mdh import completed in %s\n', my_datestr(telapsedmdh));

% if we only want the mb trace lines, swap in the trace mdh
if (~isempty(mdh_trc))
    if (readTrcRef)
        fprintf('Found MB trace reference lines --> only reading these!\n');
        mdh = mdh_trc;
    else
        fprintf('Found MB trace reference lines --> skipping these!\n');
    end
end



if (nargin >= 10)  % modification for truncation; 20170718
[mdh]=downsampled_mdh_data(mdh,[],varargin{9},varargin{10});
end


tstartmemalloc = tic;

% store k-space center line in mrprot for later use
mrprot.eja_ushKSpaceCentreLineNo = max(mdh.ushKSpaceCentreLineNo);

% now that we know how many fids we have, and how big they are,
% we can allocate memory for them
if (deOversample == true)
    OSfactor = mrprot.flReadoutOSFactor;
else
    OSfactor = 1.0;
end

% for VDxx with the channelheader concept, allocate a 3D fid array; otherwise, 2D
if (~VD11)
    nCH = 1;
else
    nCH = size(mdh.sCH_ulTypeAndChannelLength,2);
end

% find the maximum number of samples for actual image data
isACQEND     = (bitget(mdh.aulEvalInfoMask(:,1),1) == 1);
isRTFEEDBACK = ( (bitget(mdh.aulEvalInfoMask(:,1),2) == 1) & (bitget(mdh.aulEvalInfoMask(:,1),22) ~= 1) ); % not PHASCOR
isHPFEEDBACK = (bitget(mdh.aulEvalInfoMask(:,1),3) == 1);
if (nnz(isRTFEEDBACK)), fprintf('Found RTFEEDBACK lines --> skipping these!\n'); end
if (nnz(isHPFEEDBACK)), fprintf('Found HPFEEDBACK lines --> skipping these!\n'); end
maxSamplesInScan = max(mdh.ushSamplesInScan(~isACQEND & ~isRTFEEDBACK & ~isHPFEEDBACK));

idx = size(mdh.u64FidOffsets,1);

fprintf('Allocating memory for temporary raw data (%d x %d x %d complex single array)\n',maxSamplesInScan,idx,nCH);
tmpfid = complex(zeros(maxSamplesInScan,idx,nCH,'single'));
tmp = whos('tmpfid');
fprintf('Temporary raw data memory successfully allocated (%.1f MB)\n',tmp.bytes/1048576);

fprintf('Allocating memory for raw data (%d x %d x %d complex single array)\n',maxSamplesInScan/OSfactor,idx,nCH);
fid = complex(zeros(maxSamplesInScan/OSfactor,idx,nCH,'single'));
tmp = whos('fid');
fprintf('Raw data memory successfully allocated (%.1f MB)\n',tmp.bytes/1048576);

% precompute some of the interpolation parameters if necessary
lastSamples = 0; trapint = []; defaultTraj = []; destTraj = []; destTraj2 = [];
if (reInterpolate)
    % check if we really need to do this
    if (mrprot.alRegridMode(1) == 2)                     % 2 = REGRID_TRAPEZOIDAL
        % build the trapezoidal waveform
        if (interp_scheme == 2), fprintf('Reticulating splines...\n'); end
        rotrap = ones(1,mrprot.alRegridRampupTime(1)+mrprot.alRegridFlattopTime(1)+mrprot.alRegridRampdownTime(1),'single');
        roramp = single(0:1/(mrprot.alRegridRampupTime(1)+1):1);
        roramp = roramp(2:end-1);
        rotrap(1:mrprot.alRegridRampupTime(1)) = roramp(1:mrprot.alRegridRampupTime(1));
        rotrap(mrprot.alRegridRampupTime(1)+mrprot.alRegridFlattopTime(1)+1:end) = roramp(end:-1:1);
        
        % cut off the unused parts
        rotrap = rotrap(mrprot.alRegridDelaySamplesTime(1)+1:end);
        rotrap = rotrap(1:floor(mrprot.aflRegridADCDuration(1)+2)); % eja: floor added for VD11
        
        % scale to match bandwidth
        rotrap = rotrap * (length(rotrap)-1)/sum(rotrap);
        
        % integrate
        trapint = zeros(1,length(rotrap)-1,'single');
        for z=1:length(rotrap)-1
            trapint(z) = sum(rotrap(1:z));
        end
        
        % assemble the target k-space trajectory (i.e., perfect linear trajectory, time points at center of samples)
        destTraj = single(0:mrprot.aflRegridADCDuration(1)/(mrprot.alRegridDestSamples(1)+1):mrprot.aflRegridADCDuration(1));
        destTraj = destTraj(2:end-1);
        destTraj2 = destTraj(2:end-1);
        
        % precalculate the actual k-space trajectory for the first line
        nSamples = mdh.ushSamplesInScan(1);
        actualDwell = mrprot.aflRegridADCDuration(1)/double(nSamples);
        defaultTraj = interp1(1:length(trapint),trapint,1:actualDwell:mrprot.aflRegridADCDuration(1),'linear')';
        defaultTraj = defaultTraj';
        lastSamples = nSamples;
    end
end

telapsedmemalloc = toc(tstartmemalloc);
fprintf('   ** memory allocation completed in %s\n', my_datestr(telapsedmemalloc));

% call the mex to read in the FIDs
tstartdataread = tic;
rdMeas_mex_fid(outPath, tmpfid, mdh.u64FidOffsets, mdh.aulEvalInfoMask(:,1), mdh.ushSamplesInScan, mdh.ushKSpaceCentreColumn, nCH, doFlip);
telapseddataread = toc(tstartdataread);
fprintf('   ** read data completed in %s\n', my_datestr(telapseddataread));

tstartdataproc = tic;

% loop through the data and process it - only if absolutely necessary
if ( (reInterpolate) || (OSfactor > 1.0) || (doRoFT) || (doChannelScaling) )
    tremaining = -1;

    if (doChannelScaling)
        error('ERROR: this version of rdMeas does not currently support channel scaling');
    end
                
    proc_str = [];
    if ( (reInterpolate) && (mrprot.alRegridMode(1) >= 2) )
        proc_str = [proc_str '+regrid'];
    end
    if (OSfactor > 1.0)
        proc_str = [proc_str '+deOS'];
    end
    if (doRoFT)
        proc_str = [proc_str '+RoFT'];
    end
    if (~isempty(proc_str)), proc_str = [' (' proc_str ')']; end
    
    fprintf('Preprocessing image raw data%s: %10d (%3d%%)', proc_str, 0, 0);

    % steps per status update
    loop_inc = 107;
    SamplesInScan = mdh.ushSamplesInScan;
    RegridMode = mrprot.alRegridMode(1);
    RegridADCDuration = mrprot.aflRegridADCDuration(1);
    for x=1:idx % parfor should work here, but data size is often too large for available memory
        for cCH=1:nCH
            nSamples = SamplesInScan(x);
            if (nSamples)
                ctrc = tmpfid(1:nSamples,x,cCH);
                
                % interpolate the data if necessary
                if (reInterpolate)
                    % check if we really need to do this
                    if (RegridMode == 2)                     % 2 = REGRID_TRAPEZOIDAL
                        %if (findstr(dumpEvalInfoMask(mdh.aulEvalInfoMask(x,1)),'MDH_ACQEND') > 0)
                        if (isACQEND(x))
                            % skip postprocessing on this one
                        else
                            % ***** use built-in interpolation functions
                            % assemble the actual k-space trajectory for this line if it is different from the default
                            if (lastSamples ~= nSamples)
                                actualDwell = RegridADCDuration/double(nSamples);
                                actTraj = interp1(1:length(trapint),trapint,1:actualDwell:RegridADCDuration,'linear')';
                                actTraj = actTraj';
                            else
                                actTraj = defaultTraj;
                            end
                            
                            if (interp_scheme == 1)
                                % Steen modifications
                                % interpolate  Modified on 2/07/2012 to account for
                                % the difference in sampling rate from the A/D
                                % converter
                                filterf = (cumsum(destTraj2) - cumsum(actTraj'));
                                filterf = filterf / sqrt(1/length(ctrc) * sum( abs(filterf).^2) );
                                ctrc = interp1(actTraj,ctrc.*filterf',destTraj,'linear');
                            elseif (interp_scheme == 2)
                                % spline interpolation (maybe in some cases a bit better than linear)
                                ctrc = interp1(actTraj,ctrc,destTraj,'spline');
                                %elseif (interp_scheme == 3)
                                % sinc interpretation (not working well)
                                %ctrc = sinc_interp(ctrc,actTraj,destTraj);
                            else
                                % straight linear interpolation (old way)
                                ctrc = interp1(actTraj,ctrc,destTraj,'linear');
                            end
                            
                            ctrc = ctrc(2:end-1);
                            ctrc(isnan(ctrc)) = 0;
                            ctrc = ctrc';
                            ctrc = complex(real(ctrc), -imag(ctrc)); % Siemens convention
                        end
                    end
                end
                
                if (OSfactor > 1.0)
                    % de-oversample the data if desired
                    fttrc = fft(fftshift(ctrc));
                    stpos = double(nSamples)/(OSfactor*2.0);
                    fttrc = [fttrc(1:stpos); fttrc(end-stpos+1:end)];
                    if (doRoFT)
                        ctrc = fftshift(fttrc);
                    else
                        ctrc = fftshift(ifft(fttrc));
                    end
                elseif (doRoFT)
                    ctrc = fftshift(fft(fftshift(ctrc)));
                end
                
                if (nSamples/OSfactor < maxSamplesInScan)
                    ctrc_tmp = fid(:,x,cCH);
                    ctrc_tmp(1:nSamples/OSfactor) = ctrc;
                    ctrc = ctrc_tmp;
                end
                fid(:,x,cCH) = ctrc;
            end
        end
        
        % update waitbar periodically
        if (mod(x,loop_inc) == 0)
            pct_complete = x/idx;
            telapsedproc = toc(tstartdataproc);
            if (telapsedproc > 20)
                if (tremaining < 0), fprintf('                          '); end
                tremaining = telapsedproc/pct_complete - telapsedproc;
                if (tremaining > 362439), tremaining = 362439; end % max 99:99:99 for display
                tremaining_h = floor(tremaining/3600);
                tremaining_m = floor((tremaining - tremaining_h*3600)/60);
                tremaining_s = round(tremaining - tremaining_h*3600 - tremaining_m*60);
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10d (%3d%%) [est. %02d:%02d:%02d remaining]', ...
                    x, round(100*pct_complete), tremaining_h, tremaining_m, tremaining_s);
            else
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10d (%3d%%)', x, round(100*pct_complete));
            end
        end
    end
    if (tremaining < 0)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10d (%3d%%)\n', x, 100);
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10d (%3d%%)\n', x, 100);
    end

else
    fid = tmpfid;
end

telapseddataproc = toc(tstartdataproc);
fprintf('   ** process data completed in %s\n', my_datestr(telapseddataproc));

telapsed = toc(tstart);
fprintf('Success! (Elapsed time = %s)\n\n', my_datestr(telapsed));

%--------------------------------------------------------------------------
function outstr = read_cstr(fp)
% read null-terminated variable-length string from file

outstr = char(zeros(1,1000));
inchar = char(1);

idx = 1;
while (inchar ~= char(0))
    inchar = fread(fp, 1, 'uint8=>char');
    outstr(idx) = inchar;
    idx = idx + 1;
end

outstr = c_str(outstr);

%--------------------------------------------------------------------------
function outstr = my_datestr(seconds)
% convert seconds to minimal HH:MM:SS.FFF time string

if (seconds < 60)
    outstr = sprintf('%.3f s', seconds);
elseif (seconds < 3600)
    outstr = datestr(seconds/86400, 'MM:SS.FFF');
else
    outstr = datestr(seconds/86400, 'HH:MM:SS.FFF');
end