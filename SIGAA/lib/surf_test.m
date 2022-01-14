
    clear; clc;
    localdir = pwd;
    
    batchname = ['\first_batch_01\'];

locdir = strcat(localdir,batchname);

I = rgb2gray(imread(strcat(locdir, 'ptx000.bmp')));



ref_img_gray_gray = imcrop(I);
ref_pts = detectSURFFeatures(ref_img_gray_gray);
[ref_features,  ref_validPts] = extractFeatures(ref_img_gray_gray,  ref_pts);

figure; imshow(ref_img_gray_gray);
hold on; plot(ref_pts.selectStrongest(50));
%% Visual  SURF features
figure;
subplot(5,5,3); title('First Features');
for i=1:size(ref_pts)
    scale = ref_pts(i).Scale;
    image = imcrop(ref_img_gray_gray,[ref_pts(i).Location-10*scale 20*scale 20*scale]);
    subplot(5,5,i);
    imshow(image);
    hold on;
    rectangle('Position',[5*scale 5*scale 10*scale 10*scale],'Curvature',1,'EdgeColor','g');
end


%% Compare to video frame
image = imread(strcat(locdir, 'ptx100.bmp'));
I = rgb2gray(image);
% Detect features
I_pts = detectSURFFeatures(I);
[I_features, I_validPts] = extractFeatures(I, I_pts);
figure;imshow(I);
hold on; plot(I_pts.selectStrongest(50));
%% Compare card image to video frame
index_pairs = matchFeatures(ref_features, I_features);
ref_matched_pts = ref_validPts(index_pairs(:,1)).Location;
I_matched_pts = I_validPts(index_pairs(:,2)).Location;
figure, showMatchedFeatures(I, ref_img_gray_gray, I_matched_pts, ref_matched_pts, 'montage');
title('Showing all matches');



% %% Define Geometric Transformation Objects
% %gte = vision.GeometricTransformEstimator; 
% gte = estimateGeometricTransform();
% gte.Method = 'Random Sample Consensus (RANSAC)';
% [tform_matrix, inlierIdx] = step(gte, ref_matched_pts, I_matched_pts);
% ref_inlier_pts = ref_matched_pts(inlierIdx,:);
% I_inlier_pts = I_matched_pts(inlierIdx,:);
% % Draw the lines to matched points
% figure;showMatchedFeatures(image, ref_img_gray_gray, I_inlier_pts, ref_inlier_pts, 'montage');
% title('Showing match only with Inliers');
% %% Transform the corner points 
% % This will show where the object is located in the image
% tform = maketform('affine',double(tform_matrix));
% [width, height,~] = size(ref_img_gray);
% corners = [0,0;height,0;height,width;0,width];
% new_corners = tformfwd(tform, corners(:,1),corners(:,2));
% figure;imshow(image);
% patch(new_corners(:,1),new_corners(:,2),[0 1 0],'FaceAlpha',0.5);