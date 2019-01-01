% References
% 1. Class Notes
% 2. Lowe, David G. "Distinctive image features from scale-invariant keypoints." 
% International journal of computer vision 60, no. 2 (2004): 91-110.


clc;
clear all;
close all;
warning('off')

for m=1:2
    if(m==1)
        I1=rgb2gray(imresize(imread('3_1.jpg'),[1024,1024]));
        imshow(I1);
        title('Original Image1');
        I2=im2double(I1);
        [r,c]= size(I2);
    else
        I1=rgb2gray(imresize(imread('3_2.jpg'),[1024,1024]));
        I1=I1';
        imshow(I1);
        title('Original Image2');
        I2=im2double(I1);
        [r,c]= size(I2);   
    end
                                %Solution (a)
% initializing empty octave's matrix store 596
octave1=[];                                
octave2=[];
octave3=[];

%Define Variance
a=sqrt(2);
sd=1/sqrt(2);

% 1st octave
I3 = padarray(I2,[3 3],0,'both');             % zero padding
[r1,c1]= size(I3);

% Different level in First octave 
for k1=0:4                                                        
    sigma=(a^k1)*sd;
    
%generating gaussian filter coefficent
s=0;
for i=-2:2
    for j=-2:2
    K(i+3,j+3) = exp(-(i.^2+j.^2)/(2*sigma^2));
    K(i+3,j+3) = K(i+3,j+3)./(2*pi*(sigma^2));
    s=s+K(i+3,j+3); 
    end
end
K = K./s;

%convolving image with gaussian filter
I4=zeros(r,c);
s1=0;
for i = 1:r
    for j = 1:c
        sum1 = realmin;
        for k = 1:5
           for l = 1:5
               sum1=sum1+I3(i+k-1,j+l-1).*K(k,l);
               s1=s1+sum1;
           end
        end
        I4(i,j)= sum1;
    end
end
I4=I4/s1;

%show results
I5 = uint8(255 * mat2gray(I4));
figure, imshow(I5);
title('Gaussian scale-space decomposition of First Octave');
% store different scale images in octave1
octave1=[octave1 I4];
clear I4;
end

% 2nd octave
clear I3;
clear I4;
clear I5;

I3=imresize(I2,1/2);
[r,c]= size(I3);
% reduce image size by 2 in both x & y direction
I4 = padarray(I3,[3 3],0,'both');             % zero padding
[r2,c2]= size(I4);

%staring with a^2 sigma
sd1=(a.^2)*sd;

for k1=0:4                                                        
    sigma=(a^k1)*sd1;
    
    %generating gaussian filter coefficent
s=0;
for i=-2:2
    for j=-2:2
    K(i+3,j+3) = exp(-(i.^2+j.^2)/(2*sigma^2));
    K(i+3,j+3) = K(i+3,j+3)./(2*pi*(sigma^2));
    s=s+K(i+3,j+3); 
    end
end
K = K./s;
    
%convolving image with gaussian filter
I5=zeros(r,c);
s1=0;
for i = 1:r
    for j = 1:c
        sum2 = realmin;
        for k = 1:5
           for l = 1:5
               sum2=sum2+I4(i+k-1,j+l-1).*K(k,l);
               s1=s1+sum2;
           end
        end
        I5(i,j)= sum2;
    end
end
I5=I5/s1;

%show results
I6 = uint8(255 * mat2gray(I5));
figure, imshow(I6);
title('Gaussian scale-space decomposition of Second Octave');
% store different scale images in octave1
octave2=[octave2 I5];
end
% 3rd octave
 
clear I3;
clear I4;
clear I5;
clear I6;

I3=imresize(I2,1/4);
[r,c]= size(I3);
% reduce image size by 2 in both x & y direction
I4 = padarray(I3,[3 3],0,'both');             % zero padding
[r3,c3]= size(I4);
%staring with a^4 sigma
sd2=(a.^4)*sd;

for k1=0:4                                                        
    sigma=(a^k1)*sd2;

    %generating gaussian filter coefficent
    s=0;
for i=-2:2
    for j=-2:2
        K(i+3,j+3) = exp(-(i.^2+j.^2)/(2*sigma^2));
        K(i+3,j+3) = K(i+3,j+3)./(2*pi*(sigma^2));
        s=s+K(i+3,j+3); 
    end
end
K = K./s;

%convolving image with gaussian filter
I5=zeros(r,c);
s1=0;
for i = 1:r
    for j = 1:c
        sum3 = realmin;
        for k = 1:5
           for l = 1:5
              sum3=sum3+I4(i+k-1,j+l-1).*K(k,l);
              s1=s1+sum3;
           end
        end
        I5(i,j)= sum3;
    end
end
I5=I5/s1;

%show results
I6 = uint8(255 * mat2gray(I5));
figure, imshow(I6);
title('Gaussian scale-space decomposition of Third Octave');
% store different scale images in octave1
octave3=[octave3 I5];
end

%DOG - Difference of gaussian

% difference of gaussian for 1st Octave
D11=octave1(1:1024,1:1024)-octave1(1:1024,1025:2048);        
D12=octave1(1:1024,1025:2048)-octave1(1:1024,2049:3072);
D13=octave1(1:1024,2049:3072)-octave1(1:1024,3073:4096);
D14=octave1(1:1024,3073:4096)-octave1(1:1024,4097:5120);
%show results



figure,
imshow(uint8(255 * mat2gray(D11)));
title('DOG SSD of First Octave 1st Level');
figure,
imshow(uint8(255 * mat2gray(D12)));
title('DOG SSD of First Octave 2nd Level');
figure,
imshow(uint8(255 * mat2gray(D13)));
title('DOG SSD of First Octave 3rd Level');
figure,
imshow(uint8(255 * mat2gray(D14)));
title('DOG SSD of First Octave 4rd Level');

% differnce of gaussian for 2nd Octave
D21=octave2(1:512,1:512)-octave2(1:512,513:1024);
D22=octave2(1:512,513:1024)-octave2(1:512,1025:1536);
D23=octave2(1:512,1025:1536)-octave2(1:512,1537:2048);
D24=octave2(1:512,1537:2048)-octave2(1:512,2049:2560);


figure,
imshow(uint8(255 * mat2gray(D21)));
title('DOG SSD of Second Octave 1st Level');
figure,
imshow(uint8(255 * mat2gray(D22)));
title('DOG SSD of Second Octave 2nd Level');
figure,
imshow(uint8(255 * mat2gray(D23)));
title('DOG SSD of Second Octave 3rd Level');
figure,
imshow(uint8(255 * mat2gray(D24)));
title('DOG SSD of Second Octave 4rd Level');

% difference of gaussian for 3rd Octave
D31=octave3(1:256,1:256)-octave3(1:256,257:512);        
D32=octave3(1:256,257:512)-octave3(1:256,513:768);
D33=octave3(1:256,513:768)-octave3(1:256,769:1024);
D34=octave3(1:256,769:1024)-octave3(1:256,1025:1280);

figure,
imshow(uint8(255 * mat2gray(D31)));
title('DOG SSD of Third Octave 1st Level');

figure,
imshow(uint8(255 * mat2gray(D32)));
title('DOG SSD of Third Octave 2nd Level');

figure,
imshow(uint8(255 * mat2gray(D33)));
title('DOG SSD of Third Octave 3rd Level');

figure,
imshow(uint8(255 * mat2gray(D34)));
title('DOG SSD of Third Octave 4rd Level');


                                        %Solution (b)

% Key points for First Octave 
% find exterma from DOG                                      
Point00=zeros(1024,1024);
Point01=zeros(1024,1024);
Point=zeros(1024,1024);
x=0;
y=0;
z=0;


for i=2:1023
    for j=2:1023
        flag=0; % either max or min
        if (((D12(i,j)>D12(i-1,j-1))&&(D12(i,j)>D12(i-1,j)) && (D12(i,j)>D12(i-1,j+1))&&(D12(i,j)>D12(i,j-1))&&(D12(i,j)>D12(i,j+1))&&(D12(i,j)>D12(i+1,j-1)) &&(D12(i,j)>D12(i+1,j))&&(D12(i,j)>D12(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D12(i,j)>D11(i-1,j-1))&&(D12(i,j)>D11(i-1,j)) && (D12(i,j)>D11(i-1,j+1))&&(D12(i,j)>D11(i,j-1))&& D12(i,j)>D11(i,j) &&(D12(i,j)>D11(i,j+1))&&(D12(i,j)>D11(i+1,j-1)) &&(D12(i,j)>D11(i+1,j))&&(D12(i,j)>D11(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D12(i,j)>D13(i-1,j-1))&&(D12(i,j)>D13(i-1,j)) && (D12(i,j)>D13(i-1,j+1))&&(D12(i,j)>D13(i,j-1))&& D12(i,j)>D13(i,j) &&(D12(i,j)>D13(i,j+1))&&(D12(i,j)>D13(i+1,j-1)) &&(D12(i,j)>D13(i+1,j))&&(D12(i,j)>D13(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point00(i,j) = D12(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0; 
        z1=0;
        if (((D12(i,j)<D12(i-1,j-1))&&(D12(i,j)<D12(i-1,j)) && (D12(i,j)<D12(i-1,j+1))&&(D12(i,j)<D12(i,j-1))&&(D12(i,j)<D12(i,j+1))&&(D12(i,j)<D12(i+1,j-1)) &&(D12(i,j)<D12(i+1,j))&&(D12(i,j)<D12(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D12(i,j)<D11(i-1,j-1))&&(D12(i,j)<D11(i-1,j)) && (D12(i,j)<D11(i-1,j+1))&&(D12(i,j)<D11(i,j-1))&& D12(i,j)<D11(i,j) &&(D12(i,j)<D11(i,j+1))&&(D12(i,j)<D11(i+1,j-1)) &&(D12(i,j)<D11(i+1,j))&&(D12(i,j)<D11(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D12(i,j)<D13(i-1,j-1))&&(D12(i,j)<D13(i-1,j)) && (D12(i,j)<D13(i-1,j+1))&&(D12(i,j)<D13(i,j-1))&& D12(i,j)<D13(i,j) &&(D12(i,j)<D13(i,j+1))&&(D12(i,j)<D13(i+1,j-1)) &&(D12(i,j)<D13(i+1,j))&&(D12(i,j)<D13(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point00(i,j) = D12(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end


x=0;
y=0;
z=0;


for i=2:1023
    for j=2:1023
        flag=0; % either max or min
        if (((D13(i,j)>D12(i-1,j-1))&&(D13(i,j)>D12(i-1,j)) && (D13(i,j)>D12(i-1,j+1))&&(D13(i,j)>D12(i,j-1))&&(D13(i,j)>D12(i,j+1))&&(D13(i,j)>D12(i+1,j-1)) &&(D13(i,j)>D12(i+1,j))&&(D13(i,j)>D12(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D13(i,j)>D12(i-1,j-1))&&(D13(i,j)>D12(i-1,j)) && (D13(i,j)>D12(i-1,j+1))&&(D13(i,j)>D12(i,j-1))&& D13(i,j)>D12(i,j) &&(D13(i,j)>D12(i,j+1))&&(D13(i,j)>D12(i+1,j-1)) &&(D13(i,j)>D12(i+1,j))&&(D13(i,j)>D12(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D13(i,j)>D14(i-1,j-1))&&(D13(i,j)>D14(i-1,j)) && (D13(i,j)>D14(i-1,j+1))&&(D13(i,j)>D14(i,j-1))&& D13(i,j)>D14(i,j) &&(D13(i,j)>D14(i,j+1))&&(D13(i,j)>D14(i+1,j-1)) &&(D13(i,j)>D14(i+1,j))&&(D13(i,j)>D14(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point01(i,j) = D13(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0;
        z1=0;
        if (((D13(i,j)<D13(i-1,j-1))&&(D13(i,j)<D13(i-1,j)) && (D13(i,j)<D13(i-1,j+1))&&(D13(i,j)<D13(i,j-1))&&(D13(i,j)<D13(i,j+1))&&(D13(i,j)<D13(i+1,j-1)) &&(D13(i,j)<D13(i+1,j))&&(D13(i,j)<D13(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D13(i,j)<D12(i-1,j-1))&&(D13(i,j)<D12(i-1,j)) && (D13(i,j)<D12(i-1,j+1))&&(D13(i,j)<D12(i,j-1))&& D13(i,j)<D12(i,j) &&(D13(i,j)<D12(i,j+1))&&(D13(i,j)<D12(i+1,j-1)) &&(D13(i,j)<D12(i+1,j))&&(D13(i,j)<D12(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D13(i,j)<D14(i-1,j-1))&&(D13(i,j)<D14(i-1,j)) && (D13(i,j)<D14(i-1,j+1))&&(D13(i,j)<D14(i,j-1))&& D13(i,j)<D14(i,j) &&(D13(i,j)<D14(i,j+1))&&(D13(i,j)<D14(i+1,j-1)) &&(D13(i,j)<D14(i+1,j))&&(D13(i,j)<D14(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point01(i,j) = D13(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end


% Key points for Second Octave 
% find exterma from DOG                                      
Point02=zeros(512,512);
Point03=zeros(512,512);
Point1=zeros(512,512);
x=0;
y=0;
z=0;


for i=2:511
    for j=2:511
        flag=0; % either max or min
        if (((D22(i,j)>D22(i-1,j-1))&&(D22(i,j)>D22(i-1,j)) && (D22(i,j)>D22(i-1,j+1))&&(D22(i,j)>D22(i,j-1))&&(D22(i,j)>D22(i,j+1))&&(D22(i,j)>D22(i+1,j-1)) &&(D22(i,j)>D22(i+1,j))&&(D22(i,j)>D22(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D22(i,j)>D21(i-1,j-1))&&(D22(i,j)>D21(i-1,j)) && (D22(i,j)>D21(i-1,j+1))&&(D22(i,j)>D21(i,j-1))&& D22(i,j)>D21(i,j) &&(D22(i,j)>D21(i,j+1))&&(D22(i,j)>D21(i+1,j-1)) &&(D22(i,j)>D21(i+1,j))&&(D22(i,j)>D21(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D22(i,j)>D23(i-1,j-1))&&(D22(i,j)>D23(i-1,j)) && (D22(i,j)>D23(i-1,j+1))&&(D22(i,j)>D23(i,j-1))&& D22(i,j)>D23(i,j) &&(D22(i,j)>D23(i,j+1))&&(D22(i,j)>D23(i+1,j-1)) &&(D22(i,j)>D23(i+1,j))&&(D22(i,j)>D23(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point02(i,j) = D22(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0; 
        z1=0;
        if (((D22(i,j)<D22(i-1,j-1))&&(D22(i,j)<D22(i-1,j)) && (D22(i,j)<D22(i-1,j+1))&&(D22(i,j)<D22(i,j-1))&&(D22(i,j)<D22(i,j+1))&&(D22(i,j)<D22(i+1,j-1)) &&(D22(i,j)<D22(i+1,j))&&(D22(i,j)<D22(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D22(i,j)<D21(i-1,j-1))&&(D22(i,j)<D11(i-1,j)) && (D22(i,j)<D21(i-1,j+1))&&(D22(i,j)<D21(i,j-1))&& D22(i,j)<D21(i,j) &&(D22(i,j)<D21(i,j+1))&&(D22(i,j)<D21(i+1,j-1)) &&(D22(i,j)<D21(i+1,j))&&(D22(i,j)<D21(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D22(i,j)<D23(i-1,j-1))&&(D22(i,j)<D23(i-1,j)) && (D22(i,j)<D23(i-1,j+1))&&(D22(i,j)<D23(i,j-1))&& D22(i,j)<D23(i,j) &&(D22(i,j)<D23(i,j+1))&&(D22(i,j)<D23(i+1,j-1)) &&(D22(i,j)<D23(i+1,j))&&(D22(i,j)<D23(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point02(i,j) = D22(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end


x=0;
y=0;
z=0;


for i=2:511
    for j=2:511
        flag=0; % either max or min
        if (((D23(i,j)>D22(i-1,j-1))&&(D23(i,j)>D22(i-1,j)) && (D23(i,j)>D22(i-1,j+1))&&(D23(i,j)>D22(i,j-1))&&(D23(i,j)>D22(i,j+1))&&(D23(i,j)>D22(i+1,j-1)) &&(D23(i,j)>D22(i+1,j))&&(D23(i,j)>D22(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D23(i,j)>D22(i-1,j-1))&&(D23(i,j)>D22(i-1,j)) && (D23(i,j)>D22(i-1,j+1))&&(D23(i,j)>D22(i,j-1))&& D23(i,j)>D22(i,j) &&(D23(i,j)>D22(i,j+1))&&(D23(i,j)>D22(i+1,j-1)) &&(D23(i,j)>D22(i+1,j))&&(D23(i,j)>D22(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D23(i,j)>D24(i-1,j-1))&&(D23(i,j)>D24(i-1,j)) && (D23(i,j)>D24(i-1,j+1))&&(D23(i,j)>D24(i,j-1))&& D23(i,j)>D24(i,j) &&(D23(i,j)>D24(i,j+1))&&(D23(i,j)>D24(i+1,j-1)) &&(D23(i,j)>D24(i+1,j))&&(D23(i,j)>D24(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point03(i,j) = D23(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0;
        z1=0;
        if (((D23(i,j)<D23(i-1,j-1))&&(D23(i,j)<D23(i-1,j)) && (D23(i,j)<D23(i-1,j+1))&&(D23(i,j)<D23(i,j-1))&&(D23(i,j)<D23(i,j+1))&&(D23(i,j)<D23(i+1,j-1)) &&(D23(i,j)<D23(i+1,j))&&(D23(i,j)<D23(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D23(i,j)<D22(i-1,j-1))&&(D23(i,j)<D22(i-1,j)) && (D23(i,j)<D22(i-1,j+1))&&(D23(i,j)<D22(i,j-1))&& D23(i,j)<D22(i,j) &&(D23(i,j)<D22(i,j+1))&&(D23(i,j)<D22(i+1,j-1)) &&(D23(i,j)<D22(i+1,j))&&(D23(i,j)<D22(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D23(i,j)<D24(i-1,j-1))&&(D23(i,j)<D24(i-1,j)) && (D23(i,j)<D24(i-1,j+1))&&(D23(i,j)<D24(i,j-1))&& D23(i,j)<D24(i,j) &&(D23(i,j)<D24(i,j+1))&&(D23(i,j)<D24(i+1,j-1)) &&(D23(i,j)<D24(i+1,j))&&(D23(i,j)<D24(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point03(i,j) = D23(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end


% Key points for Third Octave 
% find exterma from DOG                                      
Point2=zeros(256,256);
Point04=zeros(256,256);
Point05=zeros(256,256);
x=0;
y=0;
z=0;


for i=2:255
    for j=2:255
        flag=0; % either max or min
        if (((D32(i,j)>D32(i-1,j-1))&&(D32(i,j)>D32(i-1,j)) && (D32(i,j)>D32(i-1,j+1))&&(D32(i,j)>D32(i,j-1))&&(D32(i,j)>D32(i,j+1))&&(D32(i,j)>D32(i+1,j-1)) &&(D32(i,j)>D32(i+1,j))&&(D32(i,j)>D32(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D32(i,j)>D31(i-1,j-1))&&(D32(i,j)>D31(i-1,j)) && (D32(i,j)>D31(i-1,j+1))&&(D32(i,j)>D31(i,j-1))&& D32(i,j)>D31(i,j) &&(D32(i,j)>D31(i,j+1))&&(D32(i,j)>D31(i+1,j-1)) &&(D32(i,j)>D31(i+1,j))&&(D32(i,j)>D31(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D32(i,j)>D33(i-1,j-1))&&(D32(i,j)>D33(i-1,j)) && (D32(i,j)>D33(i-1,j+1))&&(D32(i,j)>D33(i,j-1))&& D32(i,j)>D33(i,j) &&(D32(i,j)>D33(i,j+1))&&(D32(i,j)>D33(i+1,j-1)) &&(D32(i,j)>D33(i+1,j))&&(D32(i,j)>D33(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point04(i,j) = D32(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0; 
        z1=0;
        if (((D32(i,j)<D32(i-1,j-1))&&(D32(i,j)<D32(i-1,j)) && (D32(i,j)<D32(i-1,j+1))&&(D32(i,j)<D32(i,j-1))&&(D32(i,j)<D32(i,j+1))&&(D32(i,j)<D32(i+1,j-1)) &&(D32(i,j)<D32(i+1,j))&&(D32(i,j)<D32(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D32(i,j)<D31(i-1,j-1))&&(D32(i,j)<D31(i-1,j)) && (D32(i,j)<D31(i-1,j+1))&&(D32(i,j)<D31(i,j-1))&& D32(i,j)<D31(i,j) &&(D32(i,j)<D31(i,j+1))&&(D32(i,j)<D31(i+1,j-1)) &&(D32(i,j)<D31(i+1,j))&&(D32(i,j)<D31(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D32(i,j)<D33(i-1,j-1))&&(D32(i,j)<D33(i-1,j)) && (D32(i,j)<D33(i-1,j+1))&&(D32(i,j)<D33(i,j-1))&& D32(i,j)<D33(i,j) &&(D32(i,j)<D33(i,j+1))&&(D32(i,j)<D33(i+1,j-1)) &&(D32(i,j)<D33(i+1,j))&&(D32(i,j)<D33(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point04(i,j) = D32(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end


x=0;
y=0;
z=0;


for i=2:255
    for j=2:255
        flag=0; % either max or min
        if (((D33(i,j)>D32(i-1,j-1))&&(D33(i,j)>D32(i-1,j)) && (D33(i,j)>D32(i-1,j+1))&&(D33(i,j)>D32(i,j-1))&&(D33(i,j)>D32(i,j+1))&&(D33(i,j)>D32(i+1,j-1)) &&(D33(i,j)>D32(i+1,j))&&(D33(i,j)>D32(i+1,j+1))))
            x=x+1;
        else
            continue;
        end
        if x>0
            if (((D33(i,j)>D32(i-1,j-1))&&(D33(i,j)>D32(i-1,j)) && (D33(i,j)>D32(i-1,j+1))&&(D33(i,j)>D32(i,j-1))&& D33(i,j)>D32(i,j) &&(D33(i,j)>D32(i,j+1))&&(D33(i,j)>D32(i+1,j-1)) &&(D33(i,j)>D32(i+1,j))&&(D33(i,j)>D32(i+1,j+1))))
                y=y+1;
            else
                continue;
            end
        if y>0
            if (((D33(i,j)>D34(i-1,j-1))&&(D33(i,j)>D34(i-1,j)) && (D33(i,j)>D34(i-1,j+1))&&(D33(i,j)>D34(i,j-1))&& D33(i,j)>D34(i,j) &&(D33(i,j)>D34(i,j+1))&&(D33(i,j)>D34(i+1,j-1)) &&(D33(i,j)>D34(i+1,j))&&(D33(i,j)>D34(i+1,j+1))))
                z=z+1;
            else
                continue;
            end
        end
            Point05(i,j) = D33(i,j); %store Points location if it is maximum in its neighbourhood
            flag=1;
        end
    end
     
    if flag==0
        x1=0;
        y1=0;
        z1=0;
        if (((D33(i,j)<D33(i-1,j-1))&&(D33(i,j)<D33(i-1,j)) && (D33(i,j)<D33(i-1,j+1))&&(D33(i,j)<D33(i,j-1))&&(D33(i,j)<D33(i,j+1))&&(D33(i,j)<D33(i+1,j-1)) &&(D33(i,j)<D33(i+1,j))&&(D33(i,j)<D33(i+1,j+1))))
            x1=x1+1;
        else
            continue;
        end
        if x>0
            if (((D33(i,j)<D32(i-1,j-1))&&(D33(i,j)<D32(i-1,j)) && (D33(i,j)<D32(i-1,j+1))&&(D33(i,j)<D32(i,j-1))&& D33(i,j)<D32(i,j) &&(D33(i,j)<D32(i,j+1))&&(D33(i,j)<D32(i+1,j-1)) &&(D33(i,j)<D32(i+1,j))&&(D33(i,j)<D32(i+1,j+1))))
                y1=y1+1;
            else
                continue;
            end
        if y1>0
            if (((D33(i,j)<D34(i-1,j-1))&&(D33(i,j)<D34(i-1,j)) && (D33(i,j)<D34(i-1,j+1))&&(D33(i,j)<D34(i,j-1))&& D33(i,j)<D34(i,j) &&(D33(i,j)<D34(i,j+1))&&(D33(i,j)<D34(i+1,j-1)) &&(D33(i,j)<D34(i+1,j))&&(D33(i,j)<D34(i+1,j+1))))
                z1=z1+1;
            else
                continue;
            end
        end
        Point05(i,j) = D33(i,j); %store Points location if it is minimum in its neighbourhood
        end
    end
end

%Outer Rejection-  
% a)low contrast rejection and Spurious Edge Rejection
% b) Spurious Edge Detection logic
%Hessian of D

D1=[D12 D13];
D2=[D22 D23];
D3=[D32 D33];

D1= padarray(D1,[1 1],0,'both');
D2= padarray(D2,[1 1],0,'both');
D3= padarray(D3,[1 1],0,'both');

%First Octave
for i=2:1025
    for j=2:2049
        D1xx(i,j) = D1(i,j+1)-2*D1(i,j)+ D1(i,j-1); 
    end 
end

D1xx=D1xx(2:1025,2:2049);

for i=2:1025
    for j=2:2049
        D1yy(i,j) = D1(i+1,j)-2*D1(i,j)+ D1(i-1,j); 
    end
end

D1yy=D1yy(2:1025,2:2049);

for i=1:1025
    for j=2:2049
        D1x(i,j) = D1(i,j+1)- D1(i,j); 
    end 
end

D1x=D1x(2:1025,2:2049);
D1x= padarray(D1x,[1 1],0,'both');

for i=2:1025
    for j=2:2049
        D1xy(i,j) = D1x(i+1,j) - D1x(i,j); 
    end 
end

D1xy=D1xy(2:1025,2:2049);

%Second Octave

for i=2:513
    for j=2:1025
        D2xx(i,j) = D2(i,j+1)-2*D2(i,j)+ D2(i,j-1); 
    end 
end

D2xx=D2xx(2:513,2:1025);

for i=2:513
    for j=2:1025
        D2yy(i,j) = D2(i+1,j)-2*D2(i,j)+ D2(i-1,j); 
    end
end

D2yy=D2yy(2:513,2:1025);

for i=2:513
    for j=2:1025
        D2x(i,j) = D2(i,j+1)- D2(i,j); 
    end 
end

D2x=D2x(2:513,2:1025);
D2x= padarray(D2x,[1 1],0,'both');

for i=2:513
    for j=2:1025
        D2xy(i,j) = D2x(i+1,j) - D2x(i,j); 
    end 
end
D2xy=D2xy(2:513,2:1025);

%Third Octave


for i=2:257
    for j=2:513
        D3xx(i,j) = D3(i,j+1)-2*D3(i,j)+ D3(i,j-1); 
    end 
end

D3xx=D3xx(2:257,2:513);

for i=2:257
    for j=2:514
        D3yy(i,j) = D3(i+1,j)-2*D3(i,j)+ D3(i-1,j); 
    end
end

D3yy=D3yy(2:257,2:513);

for i=2:257
    for j=2:513
        D3x(i,j) = D3(i,j+1)- D3(i,j); 
    end 
end

D3x=D3x(2:257,2:513);
D3x= padarray(D3x,[1 1],0,'both');

for i=2:257
    for j=2:513
        D3xy(i,j) = D3x(i+1,j) - D3x(i,j); 
    end 
end
D3xy=D3xy(2:257,2:513);



% low contrast rejection with display Spurious edge detection

octave11 = uint8(255 * mat2gray(octave1));
octave22 = uint8(255 * mat2gray(octave2));
octave33 = uint8(255 * mat2gray(octave3));

for i=1:1024
    for j=1:1024
        tr=D1xx(i,j)+D1yy(i,j);
        det=(D1xx(i,j)*D1yy(i,j))-(D1xy(i,j).^2);
        if Point00(i,j)~=0 && abs(octave11(i,1024+j))>(0.03*255) && ((tr.^2)/det)<12.1 
            Point00(i,j)=1;
        else
            Point00(i,j)=0;
        end
        tr=D1xx(i,1024+j)+D1yy(i,1024+j);
        det=D1xx(i,1024+j)*D1yy(1024+j)-(D1xy(i,1024+j).^2);
        if Point01(i,j)~=0 && abs(octave11(i,2048+j))>(0.03*255) && ((tr.^2)/det)<12.1
            Point01(i,j)=1;
        else
            Point01(i,j)=0;
        end
    end
end

Point=Point00|Point01;

for i=1:512
    for j=1:512
        tr=D2xx(i,j)+D2yy(i,j);
        det=(D2xx(i,j)*D2yy(i,j))-(D2xy(i,j).^2);
        if Point02(i,j)~=0 && abs(octave22(i,512+j))>(0.03*255) && ((tr.^2)/det)<12.1
            Point02(i,j)=1;
        else
            Point02(i,j)=0;
        end
        tr=D2xx(i,512+j)+D2yy(i,512+j);
        det=D2xx(i,512+j)*D2yy(512+j)-(D2xy(i,512+j).^2);
        if Point03(i,j)~=0 && abs(octave22(i,1024+j))>(0.03*255) && ((tr.^2)/det)<12.1
            Point03(i,j)=1;
        else
            Point03(i,j)=0;
        end
    end
end

Point1=Point02|Point03;


for i=1:256
    for j=1:256
        tr=D3xx(i,j)+D3yy(i,j);
        det=(D3xx(i,j)*D3yy(i,j))-(D3xy(i,j).^2);
        if Point04(i,j)~=0 && abs(octave33(i,256+j))>(0.03*255) && ((tr.^2)/det)<12.1
            Point04(i,j)=1;
        else
            Point04(i,j)=0;
        end
        tr=D3xx(i,256+j)+D3yy(i,256+j);
        det=D3xx(i,256+j)*D3yy(256+j)-(D3xy(i,256+j).^2);
        if Point05(i,j)~=0 && abs(octave33(i,512+j))>(0.03*255) && ((tr.^2)/det)<12.1
            Point05(i,j)=1;
        else
            Point05(i,j)=0;
        end
    end
end

Point2=Point04|Point05;





%Display Key Points
%Key point detected on gaussian blurred images
figure,
imshow(uint8(255 * mat2gray(octave1(1:1024,1:1024))));
hold on;
for i=1:1024
    for j=1:1024
        if Point(i,j)~= 0
            plot(i, j, 'ro', 'LineWidth', 1, 'MarkerSize', 6);
        else
            continue
        end
    end
end
%title('Key point of First Octave detected on gaussian blurred images');


%Key point detected on gaussian blurred images
% figure,
% imshow(uint8(255 * mat2gray(octave2(1:512,1:512))));
% hold on;
for i=1:512
    for j=1:512
        if Point1(i,j)~= 0
            plot(2*i-1, 2*j-1, 'go', 'LineWidth', 1, 'MarkerSize', 12);
        else
            continue
        end
    end
end
%title('Key point of Second Octave detected on gaussian blurred images');


for i=1:256
    for j=1:256
        if Point2(i,j)~= 0
            plot(4*i-1, 4*j-1, 'bo', 'LineWidth', 1, 'MarkerSize', 18);
        else
            continue
        end
    end
end
title('Key point on gaussian blurred images');
                                %Solution (c)                          
%  Orientation assignment

%SmoothbyGaussianof1.5timesactulsd.
%ForeachkeypointinL,considera16x16neighbourhood ? Split into 36 bins of 10 deg each.
%Weightanglethetabymagnitudemag.
%Takeorientationgreatorthan80%ofmaxmagnitute.


I2 = padarray(I1,[8 8],0,'both');
figure,
imshow(I1);
hold on; 
%
%First Octave
for i=1:1024
    for j=1:1024
        if Point00(i,j)~=0 || Point01(i,j)~=0
            i2=padarray(I2(i:i+16,j:j+16),[2 2],0,'both');            
            img=filter1(i2,sd);
            img=img(3:20,3:20);
            mag=zeros(16,16);
            theta=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag(k-1,l-1)=sqrt(((img(k+1,l)-img(k-1,l)).^2)+((img(k,l+1)-img(k,l-1)).^2));
                    theta(k-1,l-1)=atan2((img(k+1,l)-img(k-1,l)),(img(k,l+1)-img(k,l-1)));
                end
            end
            bin=zeros(1,36);
            for u=1:16
                for v=1:16
                    if theta(u,v)>=0
                        deg=theta(u,v)*(180/pi);
                        if(deg==360)
                            deg=0;
                        end
                    else
                        deg=360+(theta(u,v)*(180/pi));
                        if(deg>351)
                            deg=0;
                        end
                    end
                    bin(floor(deg/10)+1) = bin(floor(deg/10)+1) + (mag(u,v));
                end
            end
            max1=bin(1);
            index=1;
           for o=1:8
                if(max1<bin(o))
                    max1=bin(o);
                    index=o;
                end
           end
    quiver(i,j,max1*cos((index-1)*pi/18),max1*sin((index-1)*pi/18));
           
        end
    end
end
%Second Octave
for i=1:512
    for j=1:512
        if Point02(i,j)~=0 || Point03(i,j)~=0
            i2=padarray(I2(i:i+16,j:j+16),[2 2],0,'both');            
            img=filter1(i2,a^2*sd);
            img=img(3:20,3:20);
            mag=zeros(16,16);
            theta=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag(k-1,l-1)=sqrt(((img(k+1,l)-img(k-1,l)).^2)+((img(k,l+1)-img(k,l-1)).^2));
                    theta(k-1,l-1)=atan2((img(k+1,l)-img(k-1,l)),(img(k,l+1)-img(k,l-1)));
                end
            end
            bin1=zeros(1,36);
            for u=1:16
                for v=1:16
                    if theta(u,v)>=0
                        deg=theta(u,v)*(180/pi);
                        if(deg==360)
                            deg=0;
                        end
                    else
                        deg=360+(theta(u,v)*(180/pi));
                        if(deg>351)
                           deg=0;
                        end
                    end
                    bin1(floor(deg/10)+1) = bin1(floor(deg/10)+1) + (mag(u,v));
                end
            end
            max1=bin1(1);
            index=1;
           for o=1:36
                if(max1<bin1(o))
                    max1=bin1(o);
                    index=o;
                end
           end
      
        quiver(2*i-1,2*j-1,max1*cos((index-1)*pi/18),max1*sin((index-1)*pi/18));   
           
        end
    end
end
%Third Octave
 for i=1:256
    for j=1:256 I1=I1';
        if Point04(i,j)~=0 || Point05(i,j)~=0
            i2=padarray(I2(i:i+16,j:j+16),[2 2],0,'both');            
            img=filter1(i2,a^4*sd);
            img=img(3:20,3:20);
            mag=zeros(16,16);
            theta=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag(k-1,l-1)=sqrt(((img(k+1,l)-img(k-1,l)).^2)+((img(k,l+1)-img(k,l-1)).^2));
                    theta(k-1,l-1)=atan2((img(k+1,l)-img(k-1,l)),(img(k,l+1)-img(k,l-1)));
                end
            end
            bin2=zeros(1,36);
            for u=1:16
                for v=1:16
                    if theta(u,v)>=0
                        deg=theta(u,v)*(180/pi);
                        if(deg==360)
                            deg=0;
                        end
                    else
                        deg=360+(theta(u,v)*(180/pi));
                        if(deg>351)
                               deg=0;
                        end
                    end
                    bin2(floor(deg/10)+1) = bin2(floor(deg/10)+1) + (mag(u,v));
                end
            end
            max1=bin2(1);
            index=1;
           for o=1:36
                if(max1<bin2(o))
                    max1=bin2(o);
                    index=o;
                end
           end
      
           quiver(4*i-1,4*j-1,max1*cos((index-1)*pi/18),max1*sin((index-1)*pi/18));
           
        end
    end
end

I2=I2(9:1032,9:1032);

clear i2;
                                    % Solution (d)
% Keypoint descriptor

%Convolvewith2Dgaussianofsigma=0.5x16onI(x,y) 
%Construct8binhistogramonevery4x4region.
%Concatenate16bineachofsize8.
%Descriptorsize=16*8=128dim.

I1 = padarray(I1,[8 8],0,'both');

descriptor=[];
for i=1:1024
    for j=1:1024
        if Point00(i,j)~=0 || Point01(i,j)~=0 
            i2=padarray(I1(i:i+16,j:j+16),[2 2],0,'both');
            img1=filter1(i2,0.5*16);
            img1=img1(3:20,3:20);
            mag1=zeros(16,16);
            theta1=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag1(k-1,l-1)=sqrt(((img1(k+1,l)-img1(k-1,l)).^2)+((img1(k,l+1)-img1(k,l-1)).^2));
                    theta1(k-1,l-1)=atan2((img1(k+1,l)-img1(k-1,l)),(img1(k,l+1)-img1(k,l-1)));
                end
            end
            bin3=zeros(1,8);
            for u=1:4:16
                for v=1:4:16
                    for w=0:3
                        for q=0:3
                            deg=0;
                            if theta1(u+w,v+q)>=0
                                deg=theta1(u+w,v+q)*(180/pi);
                                if(deg==360)
                                    deg=0;
                                end
                            else
                                deg=360+(theta1(u+w,v+q)*(180/pi));
                                if(deg>351)
                                    deg=0;
                                end
                            end
                            bin3(floor(deg/45)+1) = bin3(floor(deg/45)+1) + (mag1(u+w,v+q));
                            
                        end
                    end
                   descriptor=[descriptor, bin3]; 
                   bin3=zeros(1,8);
                end
            end
         descriptor=[descriptor,i,j];    
        end
    end
end

for i=1:512
    for j=1:512
        if Point01(i,j)~=0 || Point02(i,j)~=0 
            i2=padarray(I1(i:i+16,j:j+16),[2 2],0,'both');            
            img1=filter1(i2,0.5*16);
            img1=img1(3:20,3:20);
            mag1=zeros(16,16);
            theta1=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag1(k-1,l-1)=sqrt(((img1(k+1,l)-img1(k-1,l)).^2)+((img1(k,l+1)-img1(k,l-1)).^2));
                    theta1(k-1,l-1)=atan2((img1(k+1,l)-img1(k-1,l)),(img1(k,l+1)-img1(k,l-1)));
                end
            end
            bin3=zeros(1,8);
            for u=1:4:16
                for v=1:4:16
                    for w=0:3
                        for q=0:3
                            deg=0;
                            if theta1(u+w,v+q)>=0
                                deg=theta1(u+w,v+q)*(180/pi);
                                if(deg==360)
                                    deg=0;
                                end
                            else
                                deg=360+(theta1(u+w,v+q)*(180/pi));
                                if(deg>351)
                                    deg=0;
                                end
                            end
                            bin3(floor(deg/45)+1) = bin3(floor(deg/45)+1) + (mag1(u+w,v+q));
                            
                        end
                    end
                   descriptor=[descriptor, bin3]; 
                   bin3=zeros(1,8);
                end
            end
            descriptor=[descriptor,2*i-1,2*j-1];
        end
    end
end

for i=1:256
    for j=1:256
        if Point03(i,j)~=0 || Point04(i,j)~=0 
            i2=padarray(I1(i:i+16,j:j+16),[2 2],0,'both');            
            img1=filter1(i2,0.5*16);
            img1=img1(3:20,3:20);
            mag1=zeros(16,16);
            theta1=zeros(16,16);
            for k=2:17
                for l=2:17
                    mag1(k-1,l-1)=sqrt(((img1(k+1,l)-img1(k-1,l)).^2)+((img1(k,l+1)-img1(k,l-1)).^2));
                    theta1(k-1,l-1)=atan2((img1(k+1,l)-img1(k-1,l)),(img1(k,l+1)-img1(k,l-1)));
                end
            end
            bin3=zeros(1,8);
            for u=1:4:16
                for v=1:4:16
                    for w=0:3
                        for q=0:3
                            deg=0;
                            if theta1(u+w,v+q)>=0
                                deg=theta1(u+w,v+q)*(180/pi);
                                if(deg==360)
                                    deg=0;
                                end
                            else
                                deg=360+(theta1(u+w,v+q)*(180/pi));
                                if(deg>351)
                                    deg=0;
                                end
                            end
                            bin3(floor(deg/45)+1) = bin3(floor(deg/45)+1) + (mag1(u+w,v+q));
                            
                        end
                    end
                   descriptor=[descriptor, bin3]; 
                   bin3=zeros(1,8);
                end
            end
            descriptor=[descriptor,4*i-1,4*j-1];
        end
    end
end
[~,t2]=size(descriptor);
t1=t2/130;
l=1;
for i=1:t1
    for j=1:130
    des1(i,j)=descriptor(l);
    l=l+1;
    end
end
if m==1
    IMG1=I1;
    descriptor1=descriptor;
    feature1=des1';
else
    IMG2=I1;
    descriptor2=descriptor;
    feature2=des1';
end
end
%%                                  %Solution e
%Find the Euclidean distance between both feature descriptor and sort in accending
%order and plot top 20 match feature in both the images.                                   


matching=zeros(596,593);
for i=1:596
    for j=1:595 
       matching(i,j)=sqrt(sum((feature1(1:128,i)-feature2(1:128,j)).^2)); 
    end
end

[min1,col] = min(matching,[],2);

row = zeros(1,596);

for i=1:596
   row(i)=i;  
end

row=row';
find1 = [min1,row,col];

for row1 = 1:596
    for row = 1:595
        if find1(row,1) > find1(row+1,1)
           temp = find1(row,:);
           find1(row,:) = find1(row+1,:);
           find1(row+1,:) = temp;
        end
    end
end
z1=[];
z2=[];
for i=1:7:140
    x_1 = descriptor1(1,130*(find1(i,2))-1);
    y_1 = descriptor1(1,130*(find1(i,2)));
    x_2 = descriptor2(1,130*find1(i,3)-1);
    y_2 = descriptor2(1,130*find1(i,3));
    z1=[z1 [x_1 y_1]'];
    z2=[z2 [x_2 y_2]'];
    
end

figure; ax = axes;
showMatchedFeatures(IMG1(9:1032,9:1032),IMG2(9:1032,9:1032),z1',z2','montage','Parent',ax);
title(ax, 'Point matches');
legend(ax, 'Matched points 1','Matched points 2');