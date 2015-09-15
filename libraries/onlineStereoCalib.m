function [Re_cum, depthEstim, rangeMin, rangeMax, I1, I2, disparityMap, yErr] ...
    = onlineStereoCalib(I1o, I2o, Re_cum)
% Online stereocalibration algorithm for the use with FLEA cameras.
% Preliminary implementation which still needs proper handling for the
% calibration data

K_left = [3307.14974529182, 0, 705.473772999379;
   0, 3309.48543572787, 527.392642041207;
   0, 0, 1];
K_new_left = [3192.68680920959, 0, 571.05542755127;
            0, 3192.68680920959, 530.941833496094;
            0, 0, 1];
D_left = [-0.204787953655916, -0.0320519459069841, -0.000395006472233232, 0.00225328808918499, 0];

R_left = [0.999196735678056, 0.0059523359478106, 0.0396289427953972;
    -0.00597472316968856, 0.999982051494225, 0.000446511418428426;
    -0.0396255737291215, -0.000682924714446588, 0.999214365149177];

% leftP = cameraParameters('RadialDistortion', [-0.204787953655916, -0.0320519459069841, 0],...
%     'TangentialDistortion', [-0.000395006472233232, 0.00225328808918499], ...
%     'IntrinsicMatrix', leftI);

% rightP = cameraParameters('RadialDistortion', [-0.199626147736848, -0.0278720846937788, 0],...
%     'TangentialDistortion', [-0.00133984043792741, 0.00176054528915072], ...
%     'IntrinsicMatrix', rightI);

K_right = [3289.29059927731, 0, 690.312115439148;
   0, 3293.30688920301, 534.764967526671;
   0, 0, 1];

K_new_right = [3192.68680920959, 0, 571.05542755127;
            0, 3192.68680920959, 530.941833496094;
            0, 0, 1];
D_right = [-0.199626147736848, -0.0278720846937788, -0.00133984043792741, 0.00176054528915072, 0];

R_right = [0.999328717500452, 0.00208331067183515, 0.0365755956280553;
    -0.00206264782278419, 0.999997691139613, -0.000602660270738995;
    -0.0365767667086849, 0.000526813142756828, 0.999330707326184];

Base = -1758.91000893059;
BaseM = 0.556;


%%
R_left_cur = Re_cum*R_left;
R_right_cur = R_right;

%Re_cum = eye(3);

I1= uint8(opencvRectify(I1o, K_left,K_new_left, R_left_cur,D_left));
I2= uint8(opencvRectify(I2o, K_right,K_new_right, R_right_cur,D_right));



% change_current_figure(1);
% subplot(211); imagesc(I1); axis image; colormap gray; colorbar
% subplot(212); imagesc(I2); axis image; colormap gray; colorbar

% Matlabin Computer Vision System Toolboxin esimerkin mukaan:
%change_current_figure(7); imshowpair(I1,I2,'ColorChannels','red-cyan');
title('Composite Image (Red - Left Image, Cyan - Right Image)');

% Step 2: Collect Interest Points from Each Image
blobs1 = detectFASTFeatures(I1);
blobs2 = detectFASTFeatures(I2);

% Visualize the location and scale of the thirty strongest SURF features
% in I1 and I2
% change_current_figure(2);
% subplot(211); imshow(I1); hold on;
% plot(blobs1.selectStrongest(30));
% title('Thirty strongest SURF features in I1');
% 
% change_current_figure(2);
% subplot(212); imshow(I2); hold on;
% plot(blobs2.selectStrongest(30));
% title('Thirty strongest SURF features in I2');
% There may be several points with no correspondence...

% Compute SURF feature vectors (descriptors) for each blob
[features1, validBlobs1] = extractFeatures(I1, blobs1);
[features2, validBlobs2] = extractFeatures(I2, blobs2);

% Use the sum of absolute differences (SAD) metric to determine
% indices of matching features.
indexPairs = matchFeatures(features1, features2, ...
  'MatchThreshold', 5); % 103x2

% putative point correspondences
try
    matchedPoints1 = validBlobs1(indexPairs(:,1),:);
    matchedPoints2 = validBlobs2(indexPairs(:,2),:);

    % Show matching points on top of the composite image
    change_current_figure(3);
    subplot(211); showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
    legend('Matches in right image', 'Matches in left image');
    title('Original images with matched features');

    % The correctly matched points must satisfy epipolar constraints. 
    % Use the estimateFundamentalMatrix function to compute the fundamental
    % matrix and find the inliers that meet the epipolar constraint.

    [fMatrix, epipolarInliers, status] = estimateFundamentalMatrix(...
      matchedPoints1, matchedPoints2, 'Method', 'RANSAC', ...
      'NumTrials', 10000, 'DistanceThreshold', 0.01, 'Confidence', 99.99);
    %disp(fMatrix)

    inlierPoints1 = matchedPoints1(epipolarInliers, :);
    inlierPoints2 = matchedPoints2(epipolarInliers, :);
    featsL = double(inlierPoints1.Location);
    featsR = double(inlierPoints2.Location);

    initOptim = double([0, 0, 0, 0, 0]);
    UB = [0.08, 0.08, 0.08, 20, 0.01];
    LB = [-0.08, -0.08, -0.08, -20, -0.01];
    UB = [];
    LB = [];
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt');
    [minParam] = lsqnonlin(@(optim) bundleAdjust(optim, featsL, featsR), ...
        initOptim, LB, UB, OPTIONS);
    err = bundleAdjust(minParam, featsL, featsR);
    disp(['Mean error left after optimization is ', num2str(mean(err))]);

    xyDiff = mean(inlierPoints1.Location - inlierPoints2.Location);
    disp(['Difference in y dir is ', num2str(xyDiff(2)), ' pixels']);
    yErr = xyDiff(2);
    % 76 inlier points


    %change_current_figure(3); showMatchedFeatures(I1, I2, inlierPoints1, inlierPoints2);
    %legend('Inlier points in I1', 'Inlier points in I2');
    epiLines = epipolarLine(fMatrix', inlierPoints2.Location);
    points = lineToBorderPoints(epiLines, size(I1));
    
%     change_current_figure(5);
%     imagesc(I1); colormap Gray;
%     line(points(:, [1,3])', points(:, [2,4])');
%     axis image;
    % Rectify the images
%     [t1, t2] = estimateUncalibratedRectification(fMatrix, ...
%       inlierPoints1.Location, inlierPoints2.Location, size(I2));
% 
%     centroid1 = mean(inlierPoints1.Location);
%     centroid2 = mean(inlierPoints2.Location);
%     %Finding the rotation ranslation
%     cp1 = bsxfun(@minus, inlierPoints1.Location, centroid1)./3192;
%     cp2 = bsxfun(@minus, inlierPoints2.Location, centroid2)./3192;
%     H_mat = cp1'*cp2;
%     [u s v] = svd(H_mat);
%     R_p = u*v';
%     if det(R_p) < 0
%         warning('Negative determinant, reflection case')
%     end

    %Constructing the minimizing transform
    r = minParam(1);
    p = minParam(2);
    y = minParam(3);
    ty = minParam(4);
    tz = minParam(5);
    disp(['Min params: ', num2str(minParam)]);

    R = angle2dcm(r, p, y);

    T = eye(3);
    T(2:3,3) = T(2:3,3)+([ty;tz]./K_new_left(1,1));
    Re_cum = (R*T)*Re_cum;

    R_left_cur = Re_cum*R_left;
    %R_left_cur = R_left
    %R_left_cur(2,3) = Tx(3)/3192;

%             I1o = imread(fullfile(P,name_left)); % uint8
%             I1= uint8(opencvRectify(I1o, K_left,K_new_left, R_left_cur,D_left));
% 
%             name_right = ['right',name_left(5:end)];
%             I2o = imread(fullfile(P,name_right)); % uint8
%             I2= uint8(opencvRectify(I2o, K_right,K_new_right, R_right_cur,D_right));
% 
%             change_current_figure(1);
%             subplot(211); imagesc(I1); axis image; colormap gray; colorbar
%             subplot(212); imagesc(I2); axis image; colormap gray; colorbar
   % pause;

    %R_corr_cum = R_corr*R_corr_cum;
    %R_corr_cum = R_corr+R_corr_cum
    %R_left_cur = R_corr_cum.*R_left;
    %ty_cum = ty_cum +  (t_vect(2)/3192)*rGain
    %R_left_cur(2,3) = ty_cum;

%             R_left_cur(2,3) = R_left_cur(2, 3)-(xyDiff(2)/3192)/2;
%             disp(R_left_cur);
%             R_right_cur(2,3) = R_right_cur(2, 3)+(xyDiff(2)/3192)/2;

%             tform1 = projective2d(t1);
%             tform2 = projective2d(t2);
% 
%             I1Rect = imwarp(I1, tform1, 'OutputView', imref2d(size(I1)));
%             I2Rect = imwarp(I2, tform2, 'OutputView', imref2d(size(I2)));
% 
%             % transform the points to visualize them together with the
%             % rectified images
%             pts1Rect = transformPointsForward(tform1, inlierPoints1.Location);
%             pts2Rect = transformPointsForward(tform2, inlierPoints2.Location);
% 
%             %change_current_figure(4); showMatchedFeatures(I1Rect, I2Rect, pts1Rect, pts2Rect);
%             %legend('Inlier points in rectified I1', 'Inlier points in rectified I2');
% 
%             % Crop the overlapping area of the rectified images
%             Irectified = cvexTransformImagePair(I1, tform1, I2, tform2);
%             change_current_figure(3);
%             subplot(311); imshow(Irectified);
%             title('Rectified Stereo Images (Red - Left Image, Cyan - Right Image)');

    % --------------------------
    % Show points with plot3d?
    disparityMap = double(disparity(I1,I2, 'DisparityRange', [0 304], 'Method', 'SemiGlobal'));
    includeRange = 0.99;
    disparityMap(disparityMap < 0) = nan;
    depthEstim = K_new_left(1,1)*abs(BaseM) ./ disparityMap; 
    
    [rangeMin, rangeMax] = getValueRange(depthEstim, includeRange);
    
    change_current_figure(3);
    %subplot(313); imshow(disparityMap,[290 299], 'InitialMagnification', 50);
    
    subplot(212); imshow(depthEstim,[rangeMin rangeMax], 'InitialMagnification', 50);   
    colormap('jet');
    colorbar;
    title('The depth estimate [m]');

    %change_current_figure(4);
    %plot(offsets(1:loop));
    %title('Y-offsets')
catch err
    disp(err);
    warning('Something went wrong');
end



