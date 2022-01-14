
    clear; clc;
    localdir = pwd;
    
    batchname = ['\first_batch_01\'];

locdir = strcat(localdir,batchname);

J = rgb2gray(imread(strcat(locdir, 'ptx000.bmp')));

I = imcrop (J);

K = rgb2gray(imread(strcat(locdir, 'ptx100.bmp')));

imshowpair(I,K,'montage')

c = normxcorr2(I,K);
figure, surf(c), shading flat
pause;

%Find the peak in cross-correlation.

[ypeak, xpeak] = find(c==max(c(:)));
%Account for the padding that normxcorr2 adds.

yoffSet = ypeak-size(I,1);
xoffSet = xpeak-size(I,2);
%Display the matched area.

figure
imshow(K);
imrect(gca, [xoffSet+1, yoffSet+1, size(I,2), size(I,1)]);