%% virtualImageFuncDebug
function [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFuncDebug( eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters)
% Returns the virtual image coordinates for a point in eyeWorld space

options = optimset('TolFun',1e-4,'TolX',1e-4);

% Define an error function which is the distance between the nodal
% point of the camera and the point at which a ray impacts the
% plane that contains the camera, with the ray departing from the
% eyeWorld point at angle theta in the p1p2 plane. NOTE: The order
% of variables here is determined by the function that is called.
% To check:
%{
	rayTraceFuncs = compileVirtualImageFunc(createSceneGeometry());
    rayTraceFuncs.calcCameraNodeDistanceError2D.varNames
%}
errorFunc = @(theta) calcCameraNodeDistanceError2D_p1p2(...
    eyeWorldPoint, extrinsicTranslationVector, ...
    [deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion)], ...
    rotationCenters.azi([1 2]),...
    rotationCenters.ele([1 3]),...
    rotationCenters.tor([2 3]),...
    theta);
% Conduct an fminsearch to find the p1p2 theta that results in a
% ray that strikes as close as possible to the camera nodal point.
% Because the errorFunc returns nan for values very close to zero,
% we initialize the search with a slightly non-zero value (1e-4)
theta_p1p2=fminsearch(errorFunc,1e-4,options);
% Now repeat this process for in the p1p3 plane, setting the p1p2 plane to
% the theta value that was just found
errorFunc = @(theta) calcCameraNodeDistanceError3D(...
    eyeWorldPoint, extrinsicTranslationVector, ...
    [deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion)], ...
    rotationCenters.azi([1 2]),...
    rotationCenters.ele([1 3]),...
    rotationCenters.tor([2 3]),...
    theta_p1p2, theta);
% The fVal at the solution is the the total error (in mm) in both
% dimensions for intersecting the nodal point of the camera.
[theta_p1p3, nodalPointIntersectError]=fminsearch(errorFunc,1e-4,options);
% With both theta values calculated, now obtain the virtual image
% ray arising from the pupil plane that reflects the corneal optics
virtualImageRay = calcVirtualImageRay(eyeWorldPoint(1), eyeWorldPoint(2), eyeWorldPoint(3), theta_p1p2, theta_p1p3);
% Extract the origin of the ray, which is the virtual image eyeWorld point
virtualEyeWorldPoint = virtualImageRay(1,:);


end % virtualImageFunc
