function [i1,i2,F,lines] = correspondences(img1,img2)
I1 = rgb2gray(img1);
I2 = rgb2gray(img2);
points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);
[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);
indexPairs = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
% Estimate fundamental matrix between the stereo images.
[F,inliersIndex,status] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2,'Method','Norm8Point');
% figure; 
% showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
[row,col,~] = size(img2);
lamda =0.35;
 
i2 = zeros(size(img2));
I3 = applycform(img1,makecform('srgb2lab'));
I4 = applycform(img2,makecform('srgb2lab'));
lines = [];
binSize = 8 ;
magnif = 3 ;
half_patch_size = 3;    %  patch matching size

%For each point in the reference image, 
for i = half_patch_size+1:row-half_patch_size
    for j = half_patch_size+1:col-half_patch_size
        %compute epipolar line in the source image.
        templine =  epipolarLine(F',[i j]);
        points = lineToBorderPoints(templine,size(img2));
        %Draw epipolar line in source image
        %line(points(:,[1,3])',points(:,[2,4])');
        [x_coord,y_coord,token] = all_points_on_line(round(points));
        x_c = i;
        y_c = j;
        min1 = intmax('int64');
        % Search along the epipolar line in the source image to 
        %find the corresponding point of each location of the reference image.
        if(token == 0)%x y
            for l = 1:size(x_coord,2)
                if l==1
                    iMin = max(i-half_patch_size,1);
                    iMax = min(i+half_patch_size,col);
                    jMin = max(j-half_patch_size,1);
                    jMax = min(j+half_patch_size,row);
                    patch1 = I4(jMin:jMax,iMin:iMax,:);
                    l_av1 = sum(sum(patch1(:,:,1)))/(2*half_patch_size+1)^2;
                    a_av1 = sum(sum(patch1(:,:,2)))/(2*half_patch_size+1)^2;
                    b_av1 = sum(sum(patch1(:,:,3)))/(2*half_patch_size+1)^2;
                    f1 = [l_av1;a_av1;b_av1];
                    patch11 = single(rgb2gray(img2(jMin:jMax,iMin:iMax,:))) ;
                    Is1 = vl_imsmooth(patch11, sqrt((binSize/magnif)^2 - .25)) ;
                    [f1_sift1, d1] = vl_dsift(Is1, 'size', binSize) ;
                end
                iMin = max(x_coord(l)-half_patch_size,1);
                iMax = min(x_coord(l)+half_patch_size,col);
                jMin = max(y_coord(l)-half_patch_size,1);
                jMax = min(y_coord(l)+half_patch_size,row);
                patch2 = I3(iMin:iMax,jMin:jMax,:);
                l_av2 = sum(sum(patch2(:,:,1)))/(2*half_patch_size+1)^2;
                a_av2 = sum(sum(patch2(:,:,2)))/(2*half_patch_size+1)^2;
                b_av2 = sum(sum(patch2(:,:,3)))/(2*half_patch_size+1)^2;
                f2 = [l_av2;a_av2;b_av2];
                patch22 = single(rgb2gray(img1(jMin:jMax,iMin:iMax,:))) ;
                Is2 = vl_imsmooth(patch22, sqrt((binSize/magnif)^2 - .25)) ;
                [f2_sift2, d2] = vl_dsift(Is2, 'size', binSize) ;
                if sum((f1-f2).^2)+(lamda*(d1-d2.^2))<min1
                    x_c = x_coord(l);
                    y_c = y_coord(l);
                    min1 = sum((f1-f2).^2);
                end
            end
            %Create an image which is geometrically similar to the reference image using the image patches of the source image.
            if(x_c < 1)
                x_c  =1;
            elseif  x_c>row
                x_c = row;
            end
            if(y_c < 1)
                x_c  =1;
            elseif  y_c>col
                x_c = col;
            end
            i2(i,j,:) = img2(x_c,y_c,:);
        else %y   x
            for l = 1:size(x_coord,2)
                if l==1
                    iMin = max(i-half_patch_size,1);
                    iMax = min(i+half_patch_size,row);
                    jMin = max(j-half_patch_size,1);
                    jMax = min(j+half_patch_size,col);
                    patch1 = I4(jMin:jMax,iMin:iMax,:);
                    l_av1 = sum(sum(patch1(:,:,1)))/(2*half_patch_size+1)^2;
                    a_av1 = sum(sum(patch1(:,:,2)))/(2*half_patch_size+1)^2;
                    b_av1 = sum(sum(patch1(:,:,3)))/(2*half_patch_size+1)^2;
                    f1 = [l_av1;a_av1;b_av1];
                    patch11 = single(rgb2gray(img2(jMin:jMax,iMin:iMax,:))) ;
                    Is1 = vl_imsmooth(patch11, sqrt((binSize/magnif)^2 - .25)) ;
                    [f1_sift1, d1] = vl_dsift(Is1, 'size', binSize) ;
                end
                iMin = max(x_coord(l)-half_patch_size,1);
                iMax = min(x_coord(l)+half_patch_size,row);
                jMin = max(y_coord(l)-half_patch_size,1);
                jMax = min(y_coord(l)+half_patch_size,col);
                patch2 = I3(jMin:jMax,iMin:iMax,:);
                l_av2 = sum(sum(patch2(:,:,1)))/(2*half_patch_size+1)^2;
                a_av2 = sum(sum(patch2(:,:,2)))/(2*half_patch_size+1)^2;
                b_av2 = sum(sum(patch2(:,:,3)))/(2*half_patch_size+1)^2;
                f2 = [l_av2;a_av2;b_av2];
                patch22 = single(rgb2gray(img1(jMin:jMax,iMin:iMax,:))) ;
                Is2 = vl_imsmooth(patch22, sqrt((binSize/magnif)^2 - .25)) ;
                [f2_sift2, d2] = vl_dsift(Is2, 'size', binSize) ;
                if sum((f1-f2).^2)+(lamda*(d1-d2.^2))<min1
                    y_c = x_coord(l);
                    x_c = y_coord(l);
                    min1 = exp(-(sum((f1-f2).^2)));
                end
            end
            %Create an image which is geometrically similar to the reference image using the image patches of the source image.
            if(x_c < 1)
                x_c  =1;
            elseif  x_c>row
                x_c = row;
            end
            if(y_c < 1)
                x_c  =1;
            elseif  y_c>col
                x_c = col;
            end    
            i2(i,j,:) = img2(y_c,x_c,:);
        end
    end
    disp(i);
end

i1=img1;
end