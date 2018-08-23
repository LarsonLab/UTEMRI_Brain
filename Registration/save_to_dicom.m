%%
load('UTE_FinalSegs/2017-11-17_3T-DTI-volunteer_echotime_1.mat', 'movingregIC')
%%
!mkdir ./2017-11-17_3T-DTI-volunteer_echotime_1_DCM
    
info =dicominfo(['/data/larson/brain_uT2/2017-11-17_3T-DTI-volunteer/E5668/1']);
for i = 1:size(movingregIC, 3)

    info.InstanceNumber = i;
    dicomwrite(squeeze(abs(movingregIC(:,:,i,end))),['./E5668/1' num2str(i) '.DCM'] ,info);
end