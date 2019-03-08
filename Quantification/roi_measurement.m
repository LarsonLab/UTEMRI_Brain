%% ROI mean and standard deviation measurements
% obtain uT2 component values (fraction, T2, df) and DTI values for a
% given brain ROI. Requires nifti ROI files generated from FSLeyes (.nii)

% register images to T1_1mm and T1_1mm_brain in FSL

% 03/2019 Nikhil Deveshwar

addpath(genpath('/Users/nikhil/Documents/ute_t2'));


%% 
filedate = {'2018-12-10', '2019-02-06', '2018-12-10-l1','2018-06-27', '2018-06-08',...
    '2019-02-15', '2018-07-20', 'aic'}; % specify dataset
x = 7;

%% corticospinal tract (WM)roi voxel coordinates:
% 2018-12-10: (x: 54-60 y: 44-59 z: 38-63)
% 2018-12-10 right ct: (x: 31-35 y: 49-55 z: 48-54) 
    % DTI right ct: (x: 63-68 y: 106-112 z: 78-94) 
% 2019-02-06: (x: 111-116 y: 103-113 z: 75-108)
% 2019-02-06 right ct: (x: 63-68 y: 102-110 z: 90-105)


%for n = 1:length(filedate)
    % input nifti files into matlab
    ct_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_right_ct.nii']);
    ct_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_right_ct.nii']);
    ct_im_uT2_df = niftiread([filedate{x} '_uT2_df_right_ct.nii']);
    ct_im_FA = niftiread([filedate{x} '_dti_FA_right_ct.nii']);
    ct_im_MD = niftiread([filedate{x} '_dti_MD_right_ct.nii']);
    ct_im_V1 = niftiread([filedate{x} '_dti_V1_right_ct.nii']);
    ct_im_V2 = niftiread([filedate{x} '_dti_V2_right_ct.nii']);
    ct_im_V3 = niftiread([filedate{x} '_dti_V3_right_ct.nii']);

    ct_uT2_voxel_dimensions = size(ct_im_uT2_fraction);
    z_ct_uT2 = ct_uT2_voxel_dimensions(3);

    ct_DTI_voxel_dimensions = size(ct_im_FA);
    z_ct_DTI = ct_DTI_voxel_dimensions(3);

    for i = 1:z_ct_uT2
        ct_mean_uT2_fraction_array = mean(ct_im_uT2_fraction(:,:,i))/1000;
        ct_mean_uT2_T2_array = mean(ct_im_uT2_T2(:,:,i))/100;
        ct_mean_uT2_df_array = mean(ct_im_uT2_df(:,:,i))/10;
    end
    
    for j = 1:z_ct_DTI
        ct_mean_FA_array = mean(abs(ct_im_FA(:,:,j)));
        ct_mean_MD_array = mean(ct_im_MD(:,:,j));
        ct_mean_V1_array = mean(abs(ct_im_V1(:,:,j,:)));
        ct_mean_V2_array = mean(abs(ct_im_V2(:,:,j,:)));
        ct_mean_V3_array = mean(abs(ct_im_V3(:,:,j,:)));
    end


    % ct uT2 fraction
    ct_mean_uT2_fraction = mean(ct_mean_uT2_fraction_array);
    fprintf('corticospinal tract uT2 fraction mean: %f\n',ct_mean_uT2_fraction)
    ct_std_uT2_fraction = std(ct_mean_uT2_fraction_array);
    fprintf('corticospinal tract uT2 fraction std: %f\n',ct_std_uT2_fraction)

    % ct uT2 T2
    ct_mean_uT2_T2 = mean(ct_mean_uT2_T2_array);
    fprintf('corticospinal tract uT2 T2 mean: %f\n',ct_mean_uT2_T2)
    ct_std_uT2_T2 = std(ct_mean_uT2_T2_array);
    fprintf('corticospinal tract uT2 T2 std: %f\n',ct_std_uT2_T2)

    % ct uT2 df
    ct_mean_uT2_df = mean(ct_mean_uT2_df_array);
    fprintf('corticospinal tract uT2 df mean: %f\n',ct_mean_uT2_df)
    ct_std_uT2_df = std(ct_mean_uT2_df_array);
    fprintf('corticospinal tract uT2 df std: %f\n',ct_std_uT2_df)

    % FA
    ct_mean_FA = mean(ct_mean_FA_array);
    fprintf('corticospinal tract FA mean: %f\n',ct_mean_FA)
    ct_std_FA = std(ct_mean_FA_array);
    fprintf('corticospinal tract FA std: %f\n',ct_std_FA)

    % MD
    ct_mean_MD = mean(ct_mean_MD_array);
    fprintf('corticospinal tract MD mean: %f\n',ct_mean_MD)
    ct_std_MD = std(ct_mean_MD_array);
    fprintf('corticospinal tract MD std: %f\n',ct_std_MD)

    % V1
    ct_mean_V1_z = ct_mean_V1_array(:,:,1,3);
    ct_mean_V1 = mean(ct_mean_V1_z);
    fprintf('corticospinal tract V1 mean: %f\n',ct_mean_V1)
    ct_std_V1 = std(ct_mean_V1_z);
    fprintf('corticospinal tract V1 std: %f\n',ct_std_V1)

    % V2
    ct_mean_V2_x = ct_mean_V2_array(:,:,1,2);
    ct_mean_V2 = mean(ct_mean_V2_x);
    fprintf('corticospinal tract V2 mean: %f\n',ct_mean_V2)
    ct_std_V2 = std(ct_mean_V2_x);
    fprintf('corticospinal tract V2 std: %f\n',ct_std_V2)

    % V3
    ct_mean_V3_y = ct_mean_V3_array(:,:,1,3);
    ct_mean_V3 = mean(ct_mean_V3_y);
    fprintf('corticospinal tract V3 mean: %f\n',ct_mean_V3)
    ct_std_V3 = std(ct_mean_V3_y);
    fprintf('corticospinal tract V3 std: %f\n',ct_std_V3)

    ct_values = [ct_mean_uT2_fraction, ct_std_uT2_fraction; ct_mean_uT2_T2, ct_std_uT2_T2; ...
        ct_mean_uT2_df, ct_std_uT2_df; ct_mean_FA, ct_std_FA; ct_mean_MD ct_std_MD; ...
        ct_mean_V1, ct_std_V1; ct_mean_V2, ct_std_V2; ct_mean_V3, ct_std_V3];
    
    xlswrite('corticospinal tract values', double(ct_values), 'A1');

% %end

%% corpus callosum (WM)
% genu corpus callosum roi voxel coordinates: (x: 38-47 y: 74-78 z: 36-44)
% (x: 78-102 y: 147-157 z: 71-91) 

% DTI: (x:87-93 y:151-155 z:73-83)

gcc_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_gcc.nii']);
gcc_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_gcc.nii']);
gcc_im_uT2_df = niftiread([filedate{x} '_uT2_df_gcc.nii']);
gcc_im_FA = niftiread([filedate{x} '_dti_FA_gcc.nii']);
gcc_im_MD = niftiread([filedate{x} '_dti_MD_gcc.nii']);
gcc_im_V1 = niftiread([filedate{x} '_dti_V1_gcc.nii']);
gcc_im_V2 = niftiread([filedate{x} '_dti_V2_gcc.nii']);
gcc_im_V3 = niftiread([filedate{x} '_dti_V3_gcc.nii']);


gcc_uT2_voxel_dimensions = size(gcc_im_uT2_fraction);
z_gcc_uT2 = gcc_uT2_voxel_dimensions(3);

gcc_DTI_voxel_dimensions = size(gcc_im_FA);
z_gcc_DTI = gcc_DTI_voxel_dimensions(3);

for i = 1:z_gcc_uT2
    gcc_mean_uT2_fraction_array = mean(gcc_im_uT2_fraction(:,:,i))/1000;
    gcc_mean_uT2_T2_array = mean(gcc_im_uT2_T2(:,:,i))/100;
    gcc_mean_uT2_df_array = mean(gcc_im_uT2_df(:,:,i))/10;
end

for j = 1:z_gcc_DTI 
    gcc_mean_FA_array = mean(gcc_im_FA(:,:,j));
    gcc_mean_MD_array = mean(gcc_im_MD(:,:,j));
    gcc_mean_V1_array = mean(abs(gcc_im_V1(:,:,j,:)));
    gcc_mean_V2_array = mean(abs(gcc_im_V2(:,:,j,:)));
    gcc_mean_V3_array = mean(abs(gcc_im_V3(:,:,j,:)));
end

% gcc uT2 fraction
gcc_mean_uT2_fraction = mean(gcc_mean_uT2_fraction_array);
fprintf('genu corpus callosum uT2 fraction mean: %f\n',gcc_mean_uT2_fraction)
gcc_std_uT2_fraction = std(gcc_mean_uT2_fraction_array);
fprintf('genu corpus callosum uT2 fraction std: %f\n',gcc_std_uT2_fraction)

% gcc uT2 T2
gcc_mean_uT2_T2 = mean(gcc_mean_uT2_T2_array);
fprintf('genu corpus callosum uT2 T2 mean: %f\n',gcc_mean_uT2_T2)
gcc_std_uT2_T2 = std(gcc_mean_uT2_T2_array);
fprintf('genu corpus callosum uT2 T2 std: %f\n',gcc_std_uT2_T2)

% gcc uT2 df
gcc_mean_uT2_df = mean(gcc_mean_uT2_df_array);
fprintf('genu corpus callosum uT2 df mean: %f\n',gcc_mean_uT2_df)
gcc_std_uT2_df = std(gcc_mean_uT2_df_array);
fprintf('genu corpus callosum uT2 df std: %f\n',gcc_std_uT2_df)

% FA
gcc_mean_FA = mean(gcc_mean_FA_array);
fprintf('genu corpus callosum FA mean: %f\n',gcc_mean_FA)
gcc_std_FA = std(gcc_mean_FA_array);
fprintf('genu corpus callosum FA std: %f\n',gcc_std_FA)

% MD
gcc_mean_MD = mean(gcc_mean_MD_array);
fprintf('genu corpus callosum MD mean: %f\n',gcc_mean_MD)
gcc_std_MD = std(gcc_mean_MD_array);
fprintf('genu corpus callosum MD std: %f\n',gcc_std_MD)

% V1
gcc_mean_V1_y = gcc_mean_V1_array(:,:,1,1);
gcc_mean_V1 = mean(gcc_mean_V1_y);
fprintf('genu corpus callosum V1 mean: %f\n',gcc_mean_V1)
gcc_std_V1 = std(gcc_mean_V1_y);
fprintf('genu corpus callosum V1 std: %f\n',gcc_std_V1)

% V2
gcc_mean_V2_x = gcc_mean_V1_array(:,:,1,2);
gcc_mean_V2 = mean(gcc_mean_V2_x);
fprintf('genu corpus callosum V2 mean: %f\n',gcc_mean_V2)
gcc_std_V2 = std(gcc_mean_V2_x);
fprintf('genu corpus callosum V2 std: %f\n',gcc_std_V2)

% V3
gcc_mean_V3_x = gcc_mean_V1_array(:,:,1,3);
gcc_mean_V3 = mean(gcc_mean_V3_x);
fprintf('genu corpus callosum V3 mean: %f\n',gcc_mean_V3)
gcc_std_V3 = std(gcc_mean_V2_x);
fprintf('genu corpus callosum V3 std: %f\n',gcc_std_V3)

gcc_values = [gcc_mean_uT2_fraction, gcc_std_uT2_fraction; gcc_mean_uT2_T2, gcc_std_uT2_T2; ...
    gcc_mean_uT2_df, gcc_std_uT2_df; gcc_mean_FA, gcc_std_FA; gcc_mean_MD gcc_std_MD; ...
    gcc_mean_V1, gcc_std_V1; gcc_mean_V2, gcc_std_V2; gcc_mean_V3, gcc_std_V3];
 
xlswrite('gcc values', double(gcc_values), 'A1');

%% body corpus callosum (x: 42-51 y: 49-56 z: 45-48)
% (x:80-100 y: 118-126 z: 99-103)
% DTI(x: 87-94 y: 109-119 z: 97-99)

bcc_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_bcc.nii']);
bcc_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_bcc.nii']);
bcc_im_uT2_df = niftiread([filedate{x} '_uT2_df_bcc.nii']);
% bcc_im_FA = niftiread([filedate{x} '_dti_FA_bcc.nii']);
% bcc_im_MD = niftiread([filedate{x} '_dti_MD_bcc.nii']);
% bcc_im_V1 = niftiread([filedate{x} '_dti_V1_bcc.nii']);
% bcc_im_V2 = niftiread([filedate{x} '_dti_V2_bcc.nii']);
% bcc_im_V3 = niftiread([filedate{x} '_dti_V3_bcc.nii']);

bcc_uT2_voxel_dimensions = size(bcc_im_uT2_fraction);
z_bcc_uT2 = bcc_uT2_voxel_dimensions(3);

% bcc_DTI_voxel_dimensions = size(bcc_im_FA);
% z_bcc_DTI = bcc_DTI_voxel_dimensions(3);

for i = 1:z_bcc_uT2
    bcc_mean_uT2_fraction_array = mean(bcc_im_uT2_fraction(:,:,i))/1000;
    bcc_mean_uT2_T2_array = mean(bcc_im_uT2_T2(:,:,i))/100;
    bcc_mean_uT2_df_array = mean(bcc_im_uT2_df(:,:,i))/10;
end

% for j = 1:z_bcc_DTI
%     bcc_mean_FA_array = mean(abs(bcc_im_FA(:,:,j)));
%     bcc_mean_MD_array = mean(bcc_im_MD(:,:,j));
%     bcc_mean_V1_array = mean(abs(bcc_im_V1(:,:,j,:)));
%     bcc_mean_V2_array = mean(abs(bcc_im_V2(:,:,j,:)));
%     bcc_mean_V3_array = mean(abs(bcc_im_V3(:,:,j,:)));
% end

% bcc uT2 fraction
bcc_mean_uT2_fraction = mean(bcc_mean_uT2_fraction_array);
fprintf('body corpus callosum uT2 fraction mean: %f\n',bcc_mean_uT2_fraction)
bcc_std_uT2_fraction = std(bcc_mean_uT2_fraction_array);
fprintf('body corpus callosum uT2 fraction std: %f\n',bcc_std_uT2_fraction)

% bcc uT2 T2
bcc_mean_uT2_T2 = mean(bcc_mean_uT2_T2_array);
fprintf('body corpus callosum uT2 T2 mean: %f\n',bcc_mean_uT2_T2)
bcc_std_uT2_T2 = std(bcc_mean_uT2_T2_array);
fprintf('body corpus callosum uT2 T2 std: %f\n',bcc_std_uT2_T2)

% bcc uT2 df
bcc_mean_uT2_df = mean(bcc_mean_uT2_df_array);
fprintf('body corpus callosum uT2 df mean: %f\n',bcc_mean_uT2_df)
bcc_std_uT2_df = std(bcc_mean_uT2_df_array);
fprintf('body corpus callosum uT2 df std: %f\n',bcc_std_uT2_df)

% FA
bcc_mean_FA = mean(bcc_mean_FA_array);
fprintf('body corpus callosum FA mean: %f\n',bcc_mean_FA)
bcc_std_FA = std(bcc_mean_FA_array);
fprintf('body corpus callosum FA std: %f\n',bcc_std_FA)

% MD
bcc_mean_MD = mean(bcc_mean_MD_array);
fprintf('body corpus callosum MD mean: %f\n',bcc_mean_MD)
bcc_std_MD = std(bcc_mean_MD_array);
fprintf('body corpus callosum MD std: %f\n',bcc_std_MD)

% V1
bcc_mean_V1_y = bcc_mean_V1_array(:,:,1,1);
bcc_mean_V1 = mean(bcc_mean_V1_y);
fprintf('body corpus callosum V1 mean: %f\n',bcc_mean_V1)
bcc_std_V1 = std(bcc_mean_V1_y);
fprintf('body corpus callosum V1 std: %f\n',bcc_std_V1)

% V2
bcc_mean_V2_x = bcc_mean_V1_array(:,:,1,2);
bcc_mean_V2 = mean(bcc_mean_V2_x);
fprintf('body corpus callosum V2 mean: %f\n',bcc_mean_V2)
bcc_std_V2 = std(bcc_mean_V2_x);
fprintf('body corpus callosum V2 std: %f\n',bcc_std_V2)

% V3
bcc_mean_V3_x = bcc_mean_V1_array(:,:,1,3);
bcc_mean_V3 = mean(bcc_mean_V3_x);
fprintf('body corpus callosum V3 mean: %f\n',bcc_mean_V3)
bcc_std_V3 = std(bcc_mean_V2_x);
fprintf('body corpus callosum V3 std: %f\n',bcc_std_V3)

bcc_values = [bcc_mean_uT2_fraction, bcc_std_uT2_fraction; bcc_mean_uT2_T2, bcc_std_uT2_T2; ...
    bcc_mean_uT2_df, bcc_std_uT2_df; bcc_mean_FA, bcc_std_FA; bcc_mean_MD bcc_std_MD; ...
    bcc_mean_V1, bcc_std_V1; bcc_mean_V2, bcc_std_V2; bcc_mean_V3, bcc_std_V3];

xlswrite('bcc values', double(bcc_values), 'A1');

%% splenium corpus callosum (x: 42-49 y: 42-46 z: 41-48)
% (x: 80-100 y: 84-90 z: 83-89)
% DTI: (x: 76-98 y: 86-89 z: 79-87)  

scc_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_scc.nii']);
scc_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_scc.nii']);
scc_im_uT2_df = niftiread([filedate{x} '_uT2_df_scc.nii']);
scc_im_FA = niftiread([filedate{x} '_dti_FA_scc.nii']);
scc_im_MD = niftiread([filedate{x} '_dti_MD_scc.nii']);
scc_im_V1 = niftiread([filedate{x} '_dti_V1_scc.nii']);
scc_im_V2 = niftiread([filedate{x} '_dti_V2_scc.nii']);
scc_im_V3 = niftiread([filedate{x} '_dti_V3_scc.nii']);


scc_uT2_voxel_dimensions = size(scc_im_uT2_fraction);
z_scc_uT2 = scc_uT2_voxel_dimensions(3);

scc_DTI_voxel_dimensions = size(scc_im_FA);
z_scc_DTI = scc_DTI_voxel_dimensions(3);

for i = 1:z_scc_uT2
    scc_mean_uT2_fraction_array = mean(scc_im_uT2_fraction(:,:,i))/1000;
    scc_mean_uT2_T2_array = mean(scc_im_uT2_T2(:,:,i))/100;
    scc_mean_uT2_df_array = mean(scc_im_uT2_df(:,:,i))/10;
end

for j = 1:z_scc_DTI
    scc_mean_MD_array = mean(scc_im_MD(:,:,j));
    scc_mean_FA_array = mean(scc_im_FA(:,:,j));
    scc_mean_V1_array = mean(abs(scc_im_V1(:,:,j,:)));
    scc_mean_V2_array = mean(abs(scc_im_V2(:,:,j,:)));
    scc_mean_V3_array = mean(abs(scc_im_V3(:,:,j,:)));
end


% scc uT2 fraction
scc_mean_uT2_fraction = mean(scc_mean_uT2_fraction_array);
fprintf('splenium corpus callosum uT2 fraction mean: %f\n',scc_mean_uT2_fraction)
scc_std_uT2_fraction = std(scc_mean_uT2_fraction_array);
fprintf('splenium corpus callosum uT2 fraction std: %f\n',scc_std_uT2_fraction)

% scc uT2 T2
scc_mean_uT2_T2 = mean(scc_mean_uT2_T2_array);
fprintf('splenium corpus callosum uT2 T2 mean: %f\n',scc_mean_uT2_T2)
scc_std_uT2_T2 = std(scc_mean_uT2_T2_array);
fprintf('splenium corpus callosum uT2 T2 std: %f\n',scc_std_uT2_T2)

% scc uT2 df
scc_mean_uT2_df = mean(scc_mean_uT2_df_array);
fprintf('splenium corpus callosum uT2 df mean: %f\n',scc_mean_uT2_df)
scc_std_uT2_df = std(scc_mean_uT2_df_array);
fprintf('splenium corpus callosum uT2 df std: %f\n',scc_std_uT2_df)

% FA
scc_mean_FA = mean(scc_mean_FA_array);
fprintf('splenium corpus callosum FA mean: %f\n',scc_mean_FA)
scc_std_FA = std(scc_mean_FA_array);
fprintf('splenium corpus callosum FA std: %f\n',scc_std_FA)

% MD
scc_mean_MD = mean(scc_mean_MD_array);
fprintf('splenium corpus callosum MD mean: %f\n',scc_mean_MD)
scc_std_MD = std(scc_mean_MD_array);
fprintf('splenium corpus callosum MD std: %f\n',scc_std_MD)

% V1
scc_mean_V1_y = scc_mean_V1_array(:,:,1,1);
scc_mean_V1 = mean(scc_mean_V1_y);
fprintf('splenium corpus callosum V1 mean: %f\n',scc_mean_V1)
scc_std_V1 = std(scc_mean_V1_y);
fprintf('splenium corpus callosum V1 std: %f\n',scc_std_V1)

% V2
scc_mean_V2_x = scc_mean_V1_array(:,:,1,2);
scc_mean_V2 = mean(scc_mean_V2_x);
fprintf('splenium corpus callosum V2 mean: %f\n',scc_mean_V2)
scc_std_V2 = std(scc_mean_V2_x);
fprintf('splenium corpus callosum V2 std: %f\n',scc_std_V2)

% V3
scc_mean_V3_x = scc_mean_V1_array(:,:,1,3);
scc_mean_V3 = mean(scc_mean_V3_x);
fprintf('splenium corpus callosum V3 mean: %f\n',scc_mean_V3)
scc_std_V3 = std(scc_mean_V3_x);
fprintf('splenium corpus callosum V3 std: %f\n',scc_std_V3)

scc_values = [scc_mean_uT2_fraction, scc_std_uT2_fraction; scc_mean_uT2_T2, scc_std_uT2_T2; ...
    scc_mean_uT2_df, scc_std_uT2_df; scc_mean_FA, scc_std_FA; scc_mean_MD scc_std_MD; ...
    scc_mean_V1, scc_std_V1; scc_mean_V2, scc_std_V2; scc_mean_V3, scc_std_V3];

xlswrite('scc values', double(scc_values), 'A1');

%% putamen (GM) voxel coordinates: (x: 55-59 y: 61-65 z: 32-36)
% (x: 109-120 y: 130-136 z: 65-78)
% DTI: (x: 116-122 y: 119-130 z: 72-81)
pu_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_putamen.nii']);
pu_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_putamen.nii']);
pu_im_uT2_df = niftiread([filedate{x} '_uT2_df_putamen.nii']);
pu_im_FA = niftiread([filedate{x} '_dti_FA_putamen.nii']);
pu_im_MD = niftiread([filedate{x} '_dti_MD_putamen.nii']);
pu_im_V1 = niftiread([filedate{x} '_dti_V1_putamen.nii']);
pu_im_V2 = niftiread([filedate{x} '_dti_V2_putamen.nii']);
pu_im_V3 = niftiread([filedate{x} '_dti_V3_putamen.nii']);


pu_uT2_voxel_dimensions = size(pu_im_uT2_fraction);
z_pu_uT2 = pu_uT2_voxel_dimensions(3);

pu_DTI_voxel_dimensions = size(pu_im_FA);
z_pu_DTI = pu_DTI_voxel_dimensions(3);


for i = 1:z_pu_uT2
    pu_mean_uT2_fraction_array = mean(pu_im_uT2_fraction(:,:,i))/1000;
    pu_mean_uT2_T2_array = mean(pu_im_uT2_T2(:,:,i))/100;
    pu_mean_uT2_df_array = mean(pu_im_uT2_df(:,:,i))/10;
end

for j = z_pu_DTI
    pu_mean_FA_array = mean(pu_im_FA(:,:,j));
    pu_mean_MD_array = mean(pu_im_MD(:,:,j));
    pu_mean_V1_array = mean(abs(pu_im_V1(:,:,j,:)));
    pu_mean_V2_array = mean(abs(pu_im_V2(:,:,j,:)));
    pu_mean_V3_array = mean(abs(pu_im_V3(:,:,j,:)));
end

% pu uT2 fraction
pu_mean_uT2_fraction = mean(pu_mean_uT2_fraction_array);
fprintf('putamen uT2 fraction mean: %f\n',pu_mean_uT2_fraction)
pu_std_uT2_fraction = std(pu_mean_uT2_fraction_array);
fprintf('putamen uT2 fraction std: %f\n',pu_std_uT2_fraction)

% pu uT2 T2
pu_mean_uT2_T2 = mean(pu_mean_uT2_T2_array);
fprintf('putamen uT2 T2 mean: %f\n',pu_mean_uT2_T2)
pu_std_uT2_T2 = std(pu_mean_uT2_T2_array);
fprintf('putamen uT2 T2 std: %f\n',pu_std_uT2_T2)

% pu uT2 df
pu_mean_uT2_df = mean(pu_mean_uT2_df_array);
fprintf('putamen uT2 df mean: %f\n',pu_mean_uT2_df)
pu_std_uT2_df = std(pu_mean_uT2_df_array);
fprintf('putamen uT2 df std: %f\n',pu_std_uT2_df)

% FA
pu_mean_FA = mean(pu_mean_FA_array);
fprintf('putamen FA mean: %f\n',pu_mean_FA)
pu_std_FA = std(pu_mean_FA_array);
fprintf('putamen FA std: %f\n',pu_std_FA)

% MD
pu_mean_MD = mean(pu_mean_MD_array);
fprintf('putamen MD mean: %f\n',pu_mean_MD)
pu_std_MD = std(pu_mean_MD_array);
fprintf('putamen MD std: %f\n',pu_std_MD)

pu_values = [pu_mean_uT2_fraction, pu_std_uT2_fraction; pu_mean_uT2_T2, pu_std_uT2_T2; ...
    pu_mean_uT2_df, pu_std_uT2_df; pu_mean_FA, pu_std_FA; pu_mean_MD pu_std_MD;];

xlswrite('pu values', double(pu_values), 'A1');

%% caudate nucelus (GM) voxel coordinates (x: 48-51 y: 67-70 z: 38-42)
% (x: 99-106 y: 129-139 z: 80-87)


cd_im_uT2_fraction = niftiread([filedate{x} '_uT2_fraction_caudate.nii']);
cd_im_uT2_T2 = niftiread([filedate{x} '_uT2_T2_caudate.nii']);
cd_im_uT2_df = niftiread([filedate{x} '_uT2_df_caudate.nii']);
% cd_im_FA = niftiread([filedate{x} '_dti_FA_caudate.nii']);
% cd_im_MD = niftiread([filedate{x} '_dti_MD_caudate.nii']);
% cd_im_V1 = niftiread([filedate{x} '_dti_V1_caudate.nii']);
% cd_im_V2 = niftiread([filedate{x} '_dti_V2_caudate.nii']);
% cd_im_V3 = niftiread([filedate{x} '_dti_V3_caudate.nii']);


cd_uT2_voxel_dimensions = size(cd_im_uT2_fraction);
z_cd_uT2 = cd_uT2_voxel_dimensions(3);

% cd_DTI_voxel_dimensions = size(cd_im_FA);
% z_cd_DTI = cd_DTI_voxel_dimensions(3);

for i = 1:8 %z_cd_uT2
    cd_mean_uT2_fraction_array = mean(cd_im_uT2_fraction(:,:,i))/1000;
    cd_mean_uT2_T2_array = mean(cd_im_uT2_T2(:,:,i))/100;
    cd_mean_uT2_df_array = mean(cd_im_uT2_df(:,:,i))/10;
end

% for j = 1:z_cd_DTI
%     cd_mean_FA_array = mean(cd_im_FA(:,:,j));
%     cd_mean_MD_array = mean(cd_im_MD(:,:,j));
%     cd_mean_V1_array = mean(abs(cd_im_V1(:,:,j,:)));
%     cd_mean_V2_array = mean(abs(cd_im_V2(:,:,j,:)));
%     cd_mean_V3_array = mean(abs(cd_im_V3(:,:,j,:)));
% end

 

% cd uT2 fraction
cd_mean_uT2_fraction = mean(cd_mean_uT2_fraction_array);
fprintf('caudate nucleus uT2 fraction mean: %f\n',cd_mean_uT2_fraction)
cd_std_uT2_fraction = std(cd_mean_uT2_fraction_array);
fprintf('caudate nucleus uT2 fraction std: %f\n',cd_std_uT2_fraction)

% cd uT2 T2
cd_mean_uT2_T2 = mean(cd_mean_uT2_T2_array);
fprintf('caudate nucleus uT2 T2 mean: %f\n',cd_mean_uT2_T2)
cd_std_uT2_T2 = std(cd_mean_uT2_T2_array);
fprintf('caudate nucleus uT2 T2 std: %f\n',cd_std_uT2_T2)

% cd uT2 df
cd_mean_uT2_df = mean(cd_mean_uT2_df_array);
fprintf('caudate nucleus uT2 df mean: %f\n',cd_mean_uT2_df)
cd_std_uT2_df = std(cd_mean_uT2_df_array);
fprintf('caudate nucleus uT2 df std: %f\n',cd_std_uT2_df)

% FA
cd_mean_FA = mean(cd_mean_FA_array);
fprintf('caudate nucleus FA mean: %f\n',cd_mean_FA)
cd_std_FA = std(cd_mean_FA_array);
fprintf('caudate nucleus FA std: %f\n',cd_std_FA)

% MD
cd_mean_MD = mean(cd_mean_MD_array);
fprintf('caudate nucleus MD mean: %f\n',cd_mean_MD)
cd_std_MD = std(cd_mean_MD_array);
fprintf('caudate nucleus MD std: %f\n',cd_std_MD)

cd_values = [cd_mean_uT2_fraction, cd_std_uT2_fraction; cd_mean_uT2_T2, cd_std_uT2_T2; ...
    cd_mean_uT2_df, cd_std_uT2_df; cd_mean_FA, cd_std_FA; cd_mean_MD cd_std_MD;];

xlswrite('cd values', double(cd_values), 'A1');



