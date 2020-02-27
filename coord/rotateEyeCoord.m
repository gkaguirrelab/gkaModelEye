function eyeCoord = rotateEyeCoord(eyeCoord, eyePose, rotationCenters, directionFlag)
% Apply an eye rotation to an eye coordinate
%
% Syntax:
%  eyeCoord = rotateEyeCoord(eyeCoord, eyePose, rotationCenters, directionFlag)
%
% Description
%   The eye coordinate space is defined along the optical axis of the eye
%   when the eye is aligned with the optical axis of a camera. This routine
%   takes an eye coordinate and updates the value by applying an eye
%   rotation.
%
%   By default, rotation is in the forward direction specified in the
%   eyePose. If the directionFlag is set to "inverse", then the inverse
%   rotation is performed.
%
% Inputs:
%   eyeCoord              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
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
%   eyeCoord              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%

% Handle incomplete input arguments
if nargin==3
    directionFlag = 'forward';
end

% Generate the rotation matrix
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Apply the rotation, making use of the rotationCenters of the eye
if strcmp(directionFlag,'inverse')
    
    % Azimuth
    eyeCoord=eyeCoord-rotationCenters.azi;
    eyeCoord = (RotAzi'*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.azi;
    % Elevation
    eyeCoord=eyeCoord-rotationCenters.ele;
    eyeCoord = (RotEle'*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.ele;
    % Torsion
    eyeCoord=eyeCoord-rotationCenters.tor;
    eyeCoord = (RotTor'*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.tor;
    
else
    
    % Torsion
    eyeCoord=eyeCoord-rotationCenters.tor;
    eyeCoord = (RotTor*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.tor;
    % Elevation
    eyeCoord=eyeCoord-rotationCenters.ele;
    eyeCoord = (RotEle*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.ele;
    % Azimuth
    eyeCoord=eyeCoord-rotationCenters.azi;
    eyeCoord = (RotAzi*eyeCoord')';
    eyeCoord=eyeCoord+rotationCenters.azi;
    
end

end