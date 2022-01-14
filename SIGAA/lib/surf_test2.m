
    clear; clc;
    localdir = pwd;
    
    batchname = ['\first_batch_01\'];

locdir = strcat(localdir,batchname);

distorted = rgb2gray(imread(strcat(locdir, 'ptx150.bmp')));

original  = imcrop(rgb2gray(imread(strcat(locdir, 'ptx000.bmp'))));
imshow(original);
title('Base image');



figure; imshow(distorted);
title('Transformed image');


%Detect and extract features from both images.

ptsOriginal  = detectSURFFeatures(original);
ptsDistorted = detectSURFFeatures(distorted);
[featuresOriginal,validPtsOriginal] = ...
    extractFeatures(original,ptsOriginal);
[featuresDistorted,validPtsDistorted] = ...
    extractFeatures(distorted,ptsDistorted);



%Match features.

index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
figure; 
showMatchedFeatures(original,distorted,...
    matchedPtsOriginal,matchedPtsDistorted);
title('Matched SURF points,including outliers');


%Exclude the outliers, and compute the transformation matrix.

[tform,inlierPtsDistorted,inlierPtsOriginal] = ...
    estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,...
    'similarity');
figure; 

showMatchedFeatures(original,distorted,...
    inlierPtsOriginal,inlierPtsDistorted);
title('Matched inlier points');


%Recover the original image from the distorted image.

outputView = imref2d(size(original));
Ir = imwarp(distorted,tform,'OutputView',outputView);
figure; imshow(Ir); 
title('Recovered image');
 

figure
montage({original, Ir})