function ute_brain_recon(pfilename, TEfilename, offset_frequency)

if ~iscell(pfilename)
    temp = pfilename;
    clear pfilename;
    pfilename{1} = temp;
end
Npfile = length(pfilename);

TE = get_TE(TEfilename);

!mkdir tmp

for nf = 1:Npfile
imall =   precon_3dute_pfile_bartv300_allec(pfilename{nf}, 1:32, [],[], [],[],'-S -l2 -r0.01', 1, 0,0,1,1);
imallplus = precon_3dute_pfile_bartv300_allec(pfilename{nf}, 1:32, [],[],offset_frequency,[],'-S -l2 -r0.01', 1, 0,0,1,1);
save([pfilename{nf} '-csreconallec_l2_r0p01.mat'],'-v7.3', 'imall', 'imallplus','TE', 'offset_frequency');
end