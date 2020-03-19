function imagePoints = applyRadialLensDistortion(imagePointsPreDistortion,sceneGeometry)


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the sensor optical center (aka
% principal point), and the coordinates are in world units. To apply this
% distortion to our image coordinate points, we subtract the optical
% center, and then divide by fx and fy from the intrinsic matrix.

% Implement as bsxfun to avoid implicit expansion
imagePointsNormalized = bsxfun(@rdivide, ...
    bsxfun(@minus,imagePointsPreDistortion,[sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]) , ...
    [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]);

% Distortion is proportional to distance from the center of projection on
% the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(1).*radialPosition.^2 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(2).*radialPosition.^4;

imagePointsNormalizedDistorted = imagePointsNormalized.*0;
imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

% Place the distorted points back into the imagePoints vector
imagePoints = bsxfun(@plus, ...
    bsxfun(@times,imagePointsNormalizedDistorted , [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]) ,...
    [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]);

end