function [H] = homography(img1_points,img2_points)
% Estimating Homography using RANSAC
H = projective2d(eye(3));
N_pts = length(img1_points);    % total no. of points
k = 4;                          %no. of correspondences
p=0.99;                         %probability that a point is an inlier
e = 0.5;                        %outlier ratio
N_iteration = round(log10(1-p)/log10(1-(1-e)^k));

im1InlierCorrPts = [];
im2InlierCorrPts = [];
nMatchedInliers = 0; 
maxMatchedInliers = 0;

for i = 1:N_iteration
    % Random sample
    randdomIndex = randi(N_pts,k,1);
    im1pts = img1_points(randdomIndex,:);
    im2pts = img2_points(randdomIndex,:);
    % Homography Matrix estimation
    nc = length(im1pts);
    P = [im1pts(1:nc,:),ones(nc,1),zeros(nc,3),-1*im2pts(1:nc,1).*im1pts(1:nc,:),-1*im2pts(1:nc,1);
        zeros(nc,3),-1*im1pts(1:nc,:),-1*ones(nc,1),im2pts(1:nc,2).*im1pts(1:nc,:),im2pts(1:nc,2)];
    [~,S,V] = svd(P);
    sigmas = diag(S);
    if length(sigmas) >= 9
        minSigma = min(sigmas);
        [~,minSigmaCol] = find(S==minSigma);
        q = double(vpa(V(:,minSigmaCol)));
    elseif length(sigmas)<9
        q = double(vpa(V(:,9)));
    end
    A = reshape(q,[3,3])';
          
    % Geometric Error Calculation
    im1ptsFrd = [(A(1,1:2)*img1_points'+A(1,3))./(A(3,1:2)*img1_points'+A(3,3));...
    (A(2,1:2)*img1_points'+A(2,3))./(A(3,1:2)*img1_points'+A(3,3))]';
    errorForward = sum((im1ptsFrd-img2_points).^2,2).^0.5;
    totalFwdErr = sum(errorForward);
    % Backward Transformation Error
    iA=inv(A);
    im2ptsBwd = [(iA(1,1:2)*img2_points'+iA(1,3))./(iA(3,1:2)*img2_points'+iA(3,3));...
    (iA(2,1:2)*img2_points'+iA(2,3))./(iA(3,1:2)*img2_points'+iA(3,3))]';
    errorBackward = sum((im2ptsBwd-img1_points).^2,2).^0.5;
    totalBckErr = sum(errorBackward);
    % Total Geometric Error
    totalError = totalFwdErr + totalBckErr;
    sigma = sqrt(totalError/(2*N_pts));
    % Determining Threshold
    distThreshold = sqrt(5.99)*sigma;
    logic1Inlier = errorForward<distThreshold;
    logic2Inlier = errorBackward<distThreshold;
    im1Inliers = logic1Inlier.*img1_points;
    im2Inliers = logic2Inlier.*img2_points;
    matchedInliers = im1Inliers(:,1)>0 & im2Inliers(:,1)>0; 
    nMatchedInliers = nnz(matchedInliers);
    im1MatchedInliers = im1Inliers(matchedInliers,:);
    im2MatchedInliers = im2Inliers(matchedInliers,:); 
    if nMatchedInliers > maxMatchedInliers
        im1InlierCorrPts = im1MatchedInliers;
        im2InlierCorrPts = im2MatchedInliers;
        maxMatchedInliers = nMatchedInliers;
        A_ransac = A;
    end    
end
nc = length(im1InlierCorrPts);
P = [im1InlierCorrPts(1:nc,:),ones(nc,1),zeros(nc,3),...
    -1*im2InlierCorrPts(1:nc,1).*im1InlierCorrPts(1:nc,:),...
    -1*im2InlierCorrPts(1:nc,1);...
    zeros(nc,3),-1*im1InlierCorrPts(1:nc,:),-1*ones(nc,1),...
    im2InlierCorrPts(1:nc,2).*im1InlierCorrPts(1:nc,:),...
    im2InlierCorrPts(1:nc,2)];
[~,S,V] = svd(P);
sigmas = diag(S);
if length(sigmas) >= 9
    minSigma = min(sigmas);
    [~,minSigmaCol] = find(S==minSigma);
    q = double(vpa(V(:,minSigmaCol)));
    q = double(vpa(V(:,9)));
end
q = q./q(end);
H = reshape(q,[3,3])';
end