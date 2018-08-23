[X Y] = meshgrid(1:512, 1:512);
cen = [512 512]/2; rad = 75;
mask = brainmask(
img = imresize(im2double(fiximg(:,:,256)), size(mask));
figure; colormap(gray)
subplot(3,1,1)
imagesc(img); axis off
title('Image')
subplot(3,1,2)
imagesc(mask); axis off
title('Mask')
subplot(3,1,3)
imagesc(img.*mask); axis off
title('Mask Applied to Image')