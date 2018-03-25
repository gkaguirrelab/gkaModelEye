function [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFuncPreCompile( eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters)
% Calculate coord of virtual image of eyeWorldPoint
%
% Syntax:
%  virtualEyeWorldPoint = virtualImageFuncPreCompile( eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters)
%
% Description:
%   This routine calculates refraction through the optical system. In
%   application, it is never called. Instead, it is the code that is used
%   to generate a stand-alone mex file for the refraction after generating
%   stand-alone ray-tracing functions. Prior to calling this routine,
%   individual .m files must be placed on the path for each of the ray
%   tracing functions:
%   	calcCameraNodeDistanceError2D_p1p2
%       calcCameraNodeDistanceError2D_p1p3
%       calcCameraNodeDistanceError3D
%       calcVirtualImageRay
%
% SEE: compileVirtualImageFunc for details.
%

% Define an error function which is the distance between the nodal
% point of the camera and the point at which a ray impacts the
% plane that contains the camera, with the ray departing from the
% eyeWorld point at angle theta in the p1p2 plane. NOTE: The order
% of variables here is determined by the function that is called.
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
theta_p1p2=fminsearch(errorFunc,1e-4);
% Now repeat this process for a ray that varies in theta in the
% p1p3 plane
errorFunc = @(theta) calcCameraNodeDistanceError2D_p1p3(...
    eyeWorldPoint, extrinsicTranslationVector, ...
    [deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion)], ...
    rotationCenters.azi([1 2]),...
    rotationCenters.ele([1 3]),...
    rotationCenters.tor([2 3]),...
    theta);
theta_p1p3=fminsearch(errorFunc,1e-4);
% With both theta values calculated, now obtain the virtual image
% ray arising from the pupil plane that reflects the corneal optics
virtualImageRay = calcVirtualImageRay(eyeWorldPoint(1), eyeWorldPoint(2), eyeWorldPoint(3), theta_p1p2, theta_p1p3);
% Replace the original eyeWorld point with the virtual image
% eyeWorld point
virtualEyeWorldPoint = virtualImageRay(1,:);
% Calculate the total error (in mm) in both dimensions for intersecting the
% nodal point of the camera. Error values on the order of 0.01 - 0.02 are
% found across pupil points and for a range of eye rotations.
nodalPointIntersectError = ...
    calcCameraNodeDistanceError3D(...
    eyeWorldPoint, extrinsicTranslationVector, ...
    [deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion)], ...
    rotationCenters.azi([1 2]),...
    rotationCenters.ele([1 3]),...
    rotationCenters.tor([2 3]),...
    theta_p1p2, theta_p1p3);
end

