% References 
%1. Class Notes
%2. https://in.mathworks.com/matlabcentral/fileexchange/28760-2d-2d-projective-homography-3x3-estimation

%%                                Solution 1
clc;
clear all;
close all;

warning('off')

% Reading the file from folders
myfile = dir('I*'); % the folder in which your images exists
D = [myfile.isdir];
N = {myfile(D).name};

imagecell = cell(3,4);
imagecellcolor = cell(3,4);
numFiles = length(myfile);   %No.of Folder
numRows = 4;                 %No.of Files

imagecell = cell(numFiles,numRows);

for k = 1 : numFiles            
    for n = 1:numRows
    filename = sprintf('%d_%d.JPG',k,n);
    S = fullfile(N{k},filename);
    imagecellcolor{k,n} = imresize(imread(S),0.3);
    imagecell{k,n} = imresize(rgb2gray(imread(S)),0.3);
%     imageSize(n,:) = size(imagecell{k,n});
%     figure;
%     imshow(imagecell{k,n});
    end
end

% Call function Using in-bult panorama for comparision
hcellbuilt = cell(3,3);
[p1,p2,p3,hcellbuilt,tforms] = panorama(numFiles,numRows,imagecell,imagecellcolor,hcellbuilt);


%% Call function Using manual panorama for comparision
for k = 3:3            
    [panorama1,H1] = manualpanorama(imagecellcolor{k,1},imagecellcolor{k,2});
%     figure;
%     imshow(panorama1);
    [panorama2,H2] = manualpanorama(imagecellcolor{k,3},imagecellcolor{k,4});  
%     figure;
%     imshow(panorama2);
    [panorama,H3] = manualpanorama(panorama1,panorama2);
    figure;
    imshow ([imagecellcolor{k,1} panorama]);
    title('Panorama Using manual Function')
    if k==1
        panorama11=panorama;
    elseif k==2
        panorama22=panorama;
    else
        panorama33=panorama;
    end
end


%%                           Solution 2

img = cell(3,2);
imgd = cell(3,2);
h= cell(3,5);
out = cell(3,5);

h{1,4} =projective2d(eye(3));
h{2,4} =projective2d(eye(3));
h{3,4} =projective2d(eye(3));

% Reading the file from folders
img{1,1} = imread('data2/img1.jpg');
img{1,2} = imread('data2/img2.jpg');
img{2,1} = imread('data2/img3.jpg');
img{2,2} = imread('data2/img4.jpg');
img{3,1} = imread('data2/img5.jpg');
img{3,2} = imread('data2/img6.jpg');

imgd{1,1} = imread('data2/d1.jpg');
imgd{1,2} = imread('data2/d2.jpg');
imgd{2,1} = imread('data2/d3.jpg');
imgd{2,2} = imread('data2/d4.jpg');
imgd{3,1} = imread('data2/d5.jpg');
imgd{3,2} = imread('data2/d6.jpg');



for i =1:3
    p1  = detectSURFFeatures(rgb2gray(img{i,1}));
    p2 = detectSURFFeatures(rgb2gray(img{i,2}));

    [features1,  validPts1]  = extractFeatures(rgb2gray(img{i,1}),  p1);
    [features2, validPts2] = extractFeatures(rgb2gray(img{i,2}), p2);

    indexPairs = matchFeatures(features1, features2);

    matched1  = validPts1(indexPairs(:,1));
    matched2 = validPts2(indexPairs(:,2));
    
    %Quantize the depth image corresponding to the reference image into m=5 levels
    quant1 = imgd{i,1};
    quant1(imgd{i,1}<=51)                = 51;  
    quant1(imgd{i,1}<=102&imgd{i,1}>51)  = 102; 
    quant1(imgd{i,1}<=153&imgd{i,1}>102) = 153; 
    quant1(imgd{i,1}<=204&imgd{i,1}>153) = 204; 
    quant1(imgd{i,1}<=255&imgd{i,1}>154) = 255; 
    figure;
    imshow(quant1)
    
    quant2 = imgd{i,2};
    quant2(imgd{i,2}<=51)                = 51;  
    quant2(imgd{i,2}<=101&imgd{i,2}>51)  = 102; 
    quant2(imgd{i,2}<=153&imgd{i,2}>101) = 153; 
    quant2(imgd{i,2}<=204&imgd{i,2}>153) = 204; 
    quant2(imgd{i,2}<=255&imgd{i,2}>204) = 255; 
    figure;
    imshow(quant2)
    
    [r,~] = size(matched1.Location);
    p11 = [];
    p12 = [];
    p21 = [];
    p22 = [];
    p31 = [];
    p32 = [];
    p41 = [];
    p42 = [];
    p51 = [];
    p52 = [];
    
   for k =1:r
     if(quant1(round(matched1.Location(k,2)),round(matched1.Location(k,1)))==51)
         p11=[p11; round(matched1.Location(k,2)),round(matched1.Location(k,1))];
         p12=[p12; round(matched2.Location(k,2)),round(matched2.Location(k,1))];
     elseif(quant1(round(matched1.Location(k,2)),round(matched1.Location(k,1)))==102)
         p21=[p21; round(matched1.Location(k,2)),round(matched1.Location(k,1))];
         p22=[p22; round(matched2.Location(k,2)),round(matched2.Location(k,1))];
     elseif(quant1(round(matched1.Location(k,2)),round(matched1.Location(k,1)))==153)
         p31=[p31; round(matched1.Location(k,2)),round(matched1.Location(k,1))];
         p32=[p32; round(matched2.Location(k,2)),round(matched2.Location(k,1))];    
    elseif(quant1(round(matched1.Location(k,2)),round(matched1.Location(k,1)))==204)
         p41=[p41; round(matched1.Location(k,2)),round(matched1.Location(k,1))];
         p42=[p42; round(matched2.Location(k,2)),round(matched2.Location(k,1))];
     else
         p51=[p51; round(matched1.Location(k,2)),round(matched1.Location(k,1))];
         p52=[p52; round(matched2.Location(k,2)),round(matched2.Location(k,1))];
     end    
   end
   %Estimate homography matrix for each depth level of the quantized reference depth image.
   for j=1:5
      if j==1
          img1_points = p11;
          img2_points = p12;
          if(length(img1_points)~= 0)
              h{i,1}.T=homography(img1_points,img2_points);
          end
      elseif(j==2)
          img1_points = p21;
          img2_points = p22;
          if(length(img1_points)~= 0)
              h{i,2}.T=homography(img1_points,img2_points);
          end
      elseif (j==3)
          img1_points = p31;
          img2_points = p32;
          if(length(img1_points)~= 0)
              h{i,3}.T=homography(img1_points,img2_points);
          end
      elseif (j==4)
          img1_points = p41;
          img2_points = p42;
          if(length(img1_points)~= 0)
              h{i,4}.T=homography(img1_points,img2_points);
          end
      else
          img1_points = p51;
          img2_points = p52;
          if(length(img1_points)~= 0)
              h{i,5}.T=homography(img1_points,img2_points);
          end
      end
   end
   
   for j =1:5
        ih = inv(h{i,j}.T);
        h_temp = projective2d(ih');
       out{i,j} = uint8(warpping(img{i,1},h_temp));
   end
   out1 = out{i,1}+out{i,2}+out{i,3}+out{i,4};
   
   % Warp each portion of the reference RGB image corresponding to each depth level
   figure;
   imshow(out1);
   title('Warp each portion of the reference RGB image corresponding to each depth level');
   
   % Estimate a single homography matrix for the image pair and warp reference image to obtain an image similar to the source image.
   h_inbult2.T=homography(matched1.Location,matched2.Location);
   ih_inbult2.T = inv(h_inbult2.T);
   h_temp = projective2d(ih_inbult2.T');
   o = uint8(warpping(img{i,1},h_temp));
   
   figure;
   imshow(o);
   title('warp reference image to obtain an image similar to the source image') 
end