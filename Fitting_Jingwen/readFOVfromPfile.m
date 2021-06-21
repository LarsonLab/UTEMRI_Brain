function [final_FOV] = readFOVfromPfile(pfile)

header = read_MR_headers(pfile);
FOV = [header.rdb_hdr.user16, header.rdb_hdr.user17, header.rdb_hdr.user18];
undersamp = 1/header.rdb_hdr.user26;
undersamp_dc = undersamp; undersamp_imsize = undersamp;
zeropad_factor = 1;
Nramp = header.rdb_hdr.user11; %rhuser(12);
spres = header.rdb_hdr.user1; %rhuser(2);
resz_scale = header.rdb_hdr.user2; %rhuser(3);
Nprojections = header.rdb_hdr.user9; %rhuser(10);
acqs = header.rdb_hdr.user19; %rhuser(20);
shift = [header.rdb_hdr.user20, header.rdb_hdr.user21,0] / spres; %[rhuser(21:22);0].'/spres;
imsize = FOV*10/spres * zeropad_factor / undersamp_imsize;
final_imsize=round(imsize);
final_FOV = final_imsize*spres;

end