function [p1,p2,p3,hcellbuilt,tforms] = panorama(numFiles,numRows,imagecell,imagecellcolor,hcellbuilt)
%This function gives in-bult panorama for comparision for Question 1(c)
for k = 1 : numFiles   % No of folder
    points = detectSURFFeatures(imagecell{k,1});
    [features, points] = extractFeatures(imagecell{k,1}, points);
    tforms(numRows) = projective2d(eye(3));
    for n = 2:numRows        % No of files
        pointsPrevious = points;
        featuresPrevious = features;
        imageSize(n,:) = size(imagecell{k,n});
        points = detectSURFFeatures(imagecell{k,n});
        [features, points] = extractFeatures(imagecell{k,n}, points);
        indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);
        matchedPoints = points(indexPairs(:,1), :);
        matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
        tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
        tforms(n).T = tforms(n).T * tforms(n-1).T;
        hcellbuilt{k,n-1}=tforms(n).T;
    end
    % Compute the output limits  for each transform
    for i = 1:numel(tforms)
        [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
    end       
    avgXLim = mean(xlim, 2);
    [~, idx] = sort(avgXLim);
    centerIdx = floor((numel(tforms)+1)/2);
    centerImageIdx = idx(centerIdx);
    Tinv = invert(tforms(centerImageIdx));
    for i = 1:numel(tforms)
        tforms(i).T = tforms(i).T * Tinv.T;
    end
    for i = 1:numel(tforms)
       [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
    end
    maxImageSize = max(imageSize);
    % Find the minimum and maximum output limits
    xMin = min([1; xlim(:)]);
    xMax = max([maxImageSize(2); xlim(:)]);
    yMin = min([1; ylim(:)]);
    yMax = max([maxImageSize(1); ylim(:)]);
    width  = round(xMax - xMin);
    height = round(yMax - yMin);

    % Initialize the empty panorama.
    panorama = zeros([height width 3], 'like', imagecellcolor{k,n});
    blender =  vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');

    % Create a 2-D spatial reference object defining the size of the panorama.
    xLimits = [xMin xMax];
    yLimits = [yMin yMax];
    panoramaView = imref2d([height width], xLimits, yLimits);
    
    for i = 1:numRows
    warpedImage = imwarp(imagecellcolor{k,i}, tforms(i), 'OutputView', panoramaView);
    mask = imwarp(true(size(imagecellcolor{k,i},1),size(imagecellcolor{k,i},2)), tforms(i), 'OutputView', panoramaView);
    panorama = step(blender, panorama, warpedImage, mask);
    end
    if k==1
        p1=panorama;
    elseif k==2
        p2=panorama;
    else
        p3=panorama;
    end
    figure
    imshow(panorama)
    title('Panorama Using In-bult Function');
end
end