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

% Rotate each coordinate. Save the rotation matrix from the first
% coordinate to save computation time for the second coordinate
[pR, R] = rotateEyeCoord(p, eyePose, rotationCenters, directionFlag);
uR = rotateEyeCoord(u, eyePose, rotationCenters, directionFlag, R);


% Re-asemble the ray
eyeRay = quadric.normalizeRay([pR;uR-pR]')';

end