function [eyeRay, R] = rotateEyeRay(eyeRay, eyePose, rotationCenters, directionFlag, R)
% Apply an eye rotation to an eye ray
%
% Syntax:
%  [eyeRay, R] = rotateEyeRay(eyeRay, eyePose, rotationCenters, directionFlag, R)
%
% Description
%   The eye coordinate space is defined along the longitudinal axis of the
%   eye when the eye is aligned with the optical axis of a camera. This
%   routine takes an eye ray and applies an eye rotation.
%
%   By default, rotation is in the forward direction specified in the
%   eyePose. If the directionFlag is set to "inverse", then the inverse
%   rotation is performed.
%
% Inputs:
%   eyeRay                - 2x3 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The stop
%                           radius value is unused by this routine.
%   rotationCenters       - Structure. Equal to:
%                               sceneGeometry.eye.rotationCenters
%   directionFlag         - Char vector. Defaults to 'forward'.
%   R                     - Structure that defines the 3x3 rotation 
%                           matrices for each Fick angle. If passed, this
%                           rotation matrix is used instead of
%                           recalculating the matrix for the specified
%                           eyePose, saving on computation time for
%                           iterative calls to this function for the same
%                           eyePose with different eyeRay values. The
%                           field "empty" is set to true if the rotation
%                           matrices are not yet defined and consist only
%                           of nans. This convention is needed to allow
%                           code compilation.
%
% Outputs:
%   eyePoint              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   R                     - Structure that defines the 3x3 rotation 
%                           matrices for each Fick angle. 
%

% Handle incomplete input arguments
if nargin==3
    directionFlag = 'forward';
    R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
end

if nargin==4
    R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
end

if isempty(R)
    R = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);
end

% Obtain the two components of the ray as coordinates
p = eyeRay(1,:);
u = p + eyeRay(2,:);

% Rotate each coordinate. Save the rotation matrix from the first
% coordinate to save computation time for the second coordinate
[pR, R] = rotateEyeCoord(p, eyePose, rotationCenters, directionFlag, R);
uR = rotateEyeCoord(u, eyePose, rotationCenters, directionFlag, R);

% Re-asemble the ray
eyeRay = quadric.normalizeRay([pR;uR-pR]')';

end