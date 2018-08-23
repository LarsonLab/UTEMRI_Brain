%% file directories and input data



fixdir = '/data/larson/brain_uT2/2017-11-17_3T-DTI-volunteer/E5668/2';
movdir = load('/data/larson/brain_uT2/2017-11-17_3T-DTI-volunteer/P28160_20echo.7-csreconallec_l2_r0p01.mat');
echotime = 1;
id = '2017-11-17_3T-DTI-volunteer';
flip_flagfix = 1; % if 1 flip the image
flip_flagmov = 1; %1 if flip along z axis



%% load fiximg data

DLFix = dir(fullfile(fixdir,'*.DCM'));

fixfiles = cell(length(DLFix),1);
for i = 1:length(DLFix)
    fixfiles{i,1} = DLFix(i).name;
end

fixfiles = natsortfiles(fixfiles);

cd(fixdir);
fiximg = [];
for i = 1:numel(fixfiles)
    slice = dicomread(fixfiles{i});
    fiximg = cat(3, fiximg, slice);
end 

fiximg = double(fiximg);


%%
nxfix = 512; nyfix = 512; nzfix = 512;
[xfix, yfix, zfix] = ndgrid(linspace(1,size(fiximg,1),nxfix),...
                 linspace(1,size(fiximg,2),nyfix),...
                 linspace(1,size(fiximg,3),nzfix));
             
fiximg = imrotate(interp3(fiximg,xfix,yfix,zfix), -90);
if flip_flagfix
    fiximg = imrotate(fiximg,180);
end



%% load in the moving image data 
movimg = abs(movdir.imall(:,:,:,echotime));

nx = 512; ny = 512; nz = 512;
[x, y, z] = ndgrid(linspace(1,size(movimg,2),nx),...
                 linspace(1,size(movimg,1),ny),...
                 linspace(1,size(movimg,3),nz));
             
newmovimg = imrotate(interp3(movimg,x,y,z),-90);

if flip_flagmov
    newmovimg = flip(newmovimg,3);
end

 

%% perform registration

%registration only on the outline, do with binary mask


%manual adjustment of slices for each image volume 
fiximg = fiximg(:,:,100:end); %CHANGE THIS FOR EACH PATIENT
%newmovimg = newmovimg(:,:,100:end); %CHANGE THIS FOR EACH PATIENT

%create the binary mask
disp('Creating binary mask...')


fixbnry = [];
for i = 1:size(fiximg,3)
    slice = fiximg(:,:,i);
    level = mean(slice(:));
    
    slice = slice > level;
    fixbnry = cat(3,fixbnry,slice);
end

movbnry = [];
for i = 1:size(newmovimg,3)
    slice = newmovimg(:,:,i);
    level = mean(slice(:));
    
    slice = slice > level;
    movbnry = cat(3,movbnry,slice);
end


disp('Initializing optimizer and metrics...')

[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = optimizer.InitialRadius/3.5;
optimizer.MaximumIterations = 300;
optimizer.GrowthFactor = optimizer.GrowthFactor;

tform = imregtform(movbnry, fixbnry, 'rigid', optimizer, metric);
Rfixed = imref3d(size(fixbnry), 0.5195, 0.5195, 1.5);
Rmoving = imref3d(size(movbnry), .11,.9,.11);

disp('Starting Registration...')

movingreg = imwarp(movbnry, Rmoving, tform, 'OutputView', Rfixed, 'SmoothEdges', true);

movingregIC = imregister(movbnry, fixbnry, 'affine', optimizer, metric,...
    'InitialTransformation',tform);

tform_final = imregtform(movbnry, fixbnry, 'affine', optimizer, metric);


disp('Registration complete...')


%save(['/home/rcurucundhi/Documents/Rohit_Seg/UTE_FinalSegs/', ...
%    id , '_echotime_', num2str(echotime)], 'movingregIC','fiximg','tform_final');





%% code for checkerboard registration check 
image1 = double(fixbnry); 
image2 = movingregIC;
blkSize = 4;
[rows, cols, slices] = size(image1);
% blkSize = 10; % Multiple of image size

blkStart = 1:blkSize:size(image1,1);
blkEnd = blkSize:blkSize:size(image1,1);

mask1 = zeros(rows,cols);
for i = 1:2:size(blkStart,2)
    for j = 1:2:size(blkStart,2)
        mask1(blkStart(i):blkEnd(i),blkStart(j):blkEnd(j)) = 1;
    end
end

% mask2 = zeros(size(zte));
for i = 2:2:size(blkStart,2)
    for j = 2:2:size(blkStart,2)
            mask1(blkStart(i):blkEnd(i),blkStart(j):blkEnd(j)) = 1;
    end
end
mask1 = double(repmat(mask1,[1 1 slices]));
mask2 = double(~mask1);
% imagine(mask1,mask2);
regImage = mask1.*(image1) + mask2.*(image2);


save(['/home/rcurucundhi/Documents/Rohit_Seg/UTE_FinalSegs/', ...
    id , '_echotime_', num2str(echotime)], 'movingregIC','movbnry','fixbnry','tform_final'...
            ,'regImage');













