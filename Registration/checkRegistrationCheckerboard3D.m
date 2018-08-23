function regImage = checkRegistrationCheckerboard3D(image1,image2,blkSize)
% regImage = checkRegistrationCheckerboard3D(image1,image2,blkSize)
% 3d "checkerboard" to check registration
% Images need to be the same size
% Inputs
%		image1 - reference image
%		image2 - image to check with reference
%		blkSize - size of checkerboard blocks (rows and columns of image1 need to be divisible by blkSize)


% ccc;
% zte = nrrdread('ZTE_Pelvis_2mm_N4Corr.nrrd');
% ct = nrrdread('Pelvis_CT_Warped_to_ZTE.nrrd');
image1 = double(image1);
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
