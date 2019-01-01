clc;
clear all;
close all;

size = 130;

%Take a stereo image pair of static scene as the input.
I11 = imresize(imread('view0.png'),[size,size]);
I12 = imresize(imread('view1.png'),[size,size]);
I21 = imresize(imread('view2.png'),[size,size]);
I22 = imresize(imread('view3.png'),[size,size]);
I31 = imresize(imread('view4.png'),[size,size]);
I32 = imresize(imread('view5.png'),[size,size]);

[i1,i2,F1,ls1] = correspondences(I11,I12);
figure;
imshow([i1 i2 I12])
title('Image1 --- Correspondence Image --- Image2');

[i3,i4,ls2] = correspondences(I21,I22);
figure;
imshow([i3 i4 I22])
title('Image1 --- Correspondence Image --- Image2');

[i5,i6,ls3] = correspondences(I31,I32);
figure;
imshow([i5 i6 I32])
title('Image1 --- Correspondence Image --- Image2');