clc;
clear all;
close all;

%Read Image
i1=im2double(imread('butterfly.jpg'));
figure(1);
imshow(i1);
title('Original Image');

%Loop for different varience 
for i = 1:3
    if i==1
        s=1;
    elseif i==2
        s=3;
    else
        s=20;
    end
    %show results                           %Solution 1(b)
    i41=gaussian(i1,s);
    figure,
    imshow(i41);
    title(['Image Using Gaussian Filter with S.D. = ' num2str(s)]);

%     figure,
%     imshow(imgaussfilt(i1,s));
%     title(['Image Using built-in Gaussian Filter with S.D. =' num2str(s)]); 
end




%Read Image
i1=im2double(rgb2gray(imread('butterfly.jpg')));
figure(1);
imshow(i1);
title('Original Image');
[r1,c1]=size(i1);


%Kernel Size 
n1=11;   

%Variance
k=2;
s2=2;
s1=k*s2;

%Find Kernel .                          %Solution 2(a)


x = -5:5;
[Y,X] = meshgrid(x,-x);
m1 = exp( -(X.^2+Y.^2)/(2*s1^2) ); 
m1= m1./(2*pi*s1*s1);
m1 = m1 ./ sum(m1(:));


m2 = exp( -(X.^2+Y.^2)/(2*s2^2) ); 
m2= m2./(2*pi*s2*s2);
m2 = m2 ./ sum(m2(:));

DOG = m1-m2;

%Padding Image
i21 = padarray(i1,[5 5],0,'both');
[r21,c21]=size(i21);
i31=zeros(r21,c21);

%Convolution

for i = 1:r1
    for j = 1:c1
        sum2 = realmin;
        for k = 1:n1
           for l = 1:n1
               sum2=sum2+i21(i+k-1,j+l-1).*DOG(k,l);
           end
        end
        i31(i+5,j+5)= sum2;
    end
end
i31=i31/sum(i31(:));
%show results .                         %Solution 2(b)
i41 = uint8(255 * mat2gray(i31));
figure,
imshow(i41);
title('Image Using Difference of Gaussian Filter');

i31(i31 < 0) = -1;
i31(i31 > 0) = 1;
%Solution 2(c)
n=3;

i41=i31;

i51=zeros(r21,c21);

%Kernel
m3=[0 1 0;1 -4 1;0 1 0];
%Convolution
i41=double(i41);
for i = 1:r21-n
    for j = 1:c21-n
        sum3=0;
        for k = 1:n
           for l = 1:n
               sum3=sum3+i41(i+k-1,j+l-1).*m3(k,l);
           end
        end
        i51(i+1,j+1)=sum3;
    end
end

i61=zeros(r21,c21);

for i = 1:r21-1
    for j = 1:c21-1
        if (i51(i,j)* i51(i,j+1)<0)
            i61(i,j)=255;
        else
            i61(i,j)=0;
        end
        
    end
end
%show results .                         %Solution 2(c)
i71 = uint8(255 * mat2gray(i61));
figure,
imshow(i71);
title('Image Using Zero Crossings on DOG Filtered Image');



function img = gaussian(i1,var)
[r1,c1,~]=size(i1);

%Kernel Size 
n=9;
s=var;

%Find Kernel                            %Solution 1(a)
x = -4:4;
[Y,X] = meshgrid(x,-x);
global m
m = exp( -(X.^2+Y.^2)/(2*s^2) ); 
m= m./sqrt((2*pi*s^2));
m = m ./ sum(m(:));

%Padding Image
i21 = padarray(i1(:,:,1),[4 4],0,'both');
i22 = padarray(i1(:,:,2),[4 4],0,'both');
i23 = padarray(i1(:,:,3),[4 4],0,'both');

[r2,c2, t1]=size(i21);
i31=zeros(r2,c2);
i32=zeros(r2,c2);
i33=zeros(r2,c2);

% show red,blue, green components of image
% a = zeros(size(i1, 1), size(i1, 2));
% just_red = cat(3, i1(:,:,1), a, a);
% just_green = cat(3, a, i1(:,:,2), a);
% just_blue = cat(3, a, a, i1(:,:,3));


%Convolution

for i = 1:r1
    for j = 1:c1
        sum11=0;
        sum12=0;
        sum13=0;
        for k = 1:n
           for l = 1:n
               sum11=sum11+i21(i+k-1,j+l-1).*m(k,l);
               sum12=sum12+i22(i+k-1,j+l-1).*m(k,l);
               sum13=sum13+i23(i+k-1,j+l-1).*m(k,l);
           end
        end
        i31(i+4,j+4)=sum11/(n*n);
        i32(i+4,j+4)=sum12/(n*n);
        i33(i+4,j+4)=sum13/(n*n);
    end
end

%show results                           %Solution 1(b)
rgbImage = cat(3, i31, i32, i33);
img = uint8(255*mat2gray(rgbImage));
end

