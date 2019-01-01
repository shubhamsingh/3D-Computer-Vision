function [panorama,H] = manualpanorama(I1,I2)
img1 = im2double(rgb2gray(I1));
img2 = im2double(rgb2gray(I2));
points1 = detectSURFFeatures(img1);
points2 = detectSURFFeatures(img2);
features1 = extractFeatures(img1,points1);
features2 = extractFeatures(img2,points2);
indexPairs = matchFeatures(features1,features2,'Unique',true);
matchedPoints1 = points1(indexPairs(:,1));
matchedPoints2 = points2(indexPairs(:,2));
img1_points = matchedPoints1.Location;
img2_points = matchedPoints2.Location;


% Estimating Homography using RANSAC
H = homography(img1_points,img2_points);
iH = inv(H);
H = projective2d(iH');

% manual warpping
im2_transformed = warpping(I2,H);

im1_size = size(I1);
im2_transformed_size = size(im2_transformed);
padsize = im2_transformed_size-im1_size;
img1_expanded = padarray(I1,abs(padsize),'post');

%manual Blending
x_overlap=[0.121;10.23];
overlapleft = round(x_overlap(1));
overlapright = round(x_overlap(2));
diff = overlapright-overlapleft;
ramp = [zeros(1,overlapleft-1),0:1/diff:1,ones(1,im2_transformed_size(2)-overlapright)];
[~,c,~]= size(im2_transformed);
im2_blend = im2_transformed .* repmat(ramp(1,1:c),im2_transformed_size(1),1);
ramp_for_im1 = 1 - ramp;
[~,c,~]= size(img1_expanded);
im1_blend = img1_expanded .* uint8(repmat(ramp_for_im1(1,1:c),im2_transformed_size(1),1));
panorama = im1_blend+uint8(im2_blend);
panorama = I1+panorama;
% figure
% imshow(panorama)
end