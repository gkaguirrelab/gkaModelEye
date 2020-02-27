function eyeRay = rotateEyeRay(eyeRay, eyePose, rotationCenters, directionFlag)
% Apply an eye rotation to an eye ray
%
% Syntax:
%  eyeRay = rotateEyeRay(eyeRay, eyePose, rotationCenters, directionFlag)
%
% Description
%   The eye coordinate space is defined along the optical axis of the eye
%   when the eye is aligned with the optical axis of a camera. This routine
%   takes an eye ray and applies an eye rotation.
%
%   By default, rotation is in the forward direction specified in the
%   eyePose. If the directionFlag is set to "inverse", then the inverse
%   rotation is performed.
%
% Inputs:
%   eyeRay                - 2x3 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%                           dimensions p1, p2, p3.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The stop
%                           radius value is unused by this routine.
%   rotationCenters       - Structure. Equal to:
%                               sceneGeometry.eye.rotationCenters
%   directionFlag         - Char vector. Defaults to 'forward'.
%
% Outputs:
%   eyePoint              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%

% Handle incomplete input arguments
if nargin==3
    directionFlag = 'forward';
end

% Obtain the two components of the ray as coordinates
p = eyeRay(1,:);
u = p + eyeRay(2,:);

% Generate the rotation matrix
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Apply the rotation (or inverse rotation if requested), making use of the
% rotationCenters of the eye

if strcmp(directionFlag,'inverse')
    
    % Azimuth
    p=p-rotationCenters.azi;
    p = (RotAzi'*p')';
    p=p+rotationCenters.azi;
    % Elevation
    p=p-rotationCenters.ele;
    p = (RotEle'*p')';
    p=p+rotationCenters.ele;
    % Torsion
    p=p-rotationCenters.tor;
    p = (RotTor'*p')';
    p=p+rotationCenters.tor;
    
    % Azimuth
    u=u-rotationCenters.azi;
    u = (RotAzi'*u')';
    u=u+rotationCenters.azi;
    % Elevation
    u=u-rotationCenters.ele;
    u = (RotEle'*u')';
    u=u+rotationCenters.ele;
    % Torsion
    u=u-rotationCenters.tor;
    u = (RotTor'*u')';
    u=u+rotationCenters.tor;
    
else
    
    % Torsion
    p=p-rotationCenters.tor;
    p = (RotTor*p')';
    p=p+rotationCenters.tor;
    % Elevation
    p=p-rotationCenters.ele;
    p = (RotEle*p')';
    p=p+rotationCenters.ele;
    % Azimuth
    p=p-rotationCenters.azi;
    p = (RotAzi*p')';
    p=p+rotationCenters.azi;
    
    % Torsion
    u=u-rotationCenters.tor;
    u = (RotTor*u')';
    u=u+rotationCenters.tor;
    % Elevation
    u=u-rotationCenters.ele;
    u = (RotEle*u')';
    u=u+rotationCenters.ele;
    % Azimuth
    u=u-rotationCenters.azi;
    u = (RotAzi*u')';
    u=u+rotationCenters.azi;
end

% Re-asemble the ray
eyeRay = quadric.normalizeRay([p;u-p]')';

end