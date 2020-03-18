function [eyeCoord, R] = rotateEyeCoord(eyeCoord, eyePose, rotationCenters, directionFlag, R)
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
%   The eye rotation is defined in the "Fick" coordinate frame. In this
%   coordinate frame, the eye undergoes azimuthal rotation around a fixed
%   vertical axies. There is then an elevational rotation around the
%   rotated horizontal axis, and finally a torsional rotation around the
%   rotated optical axis of the eye. We implement here these rotations
%   around a head-fixed, non-rotating set of axes. Therefore, the order of
%   rotations is reversed (tor-ele-azi), but this is still the Fick
%   coordinate system.
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
% Examples:
%{
    % Test that forward and inverse rotation recovers the same value
    eye = modelEyeParameters();
    eyeCoordSource = [0, 5, -2.5];
    eyePose = [-10, 15, 3];
    ecR = rotateEyeCoord(eyeCoordSource, eyePose, eye.rotationCenters,'forward');
    eyeCoordRecovered = rotateEyeCoord(ecR, eyePose, eye.rotationCenters,'inverse');
    assert(max(abs((eyeCoordSource - eyeCoordRecovered))) < 1e-12)
%}

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


%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to rightward, upward and clockwise (as seen
% from the perspective of the subject) eye movements
if R.empty
    R.azi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
    R.ele = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
    R.tor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];
    R.empty = false;
end

% We shift the points to each rotation center, rotate, shift back, and
% repeat. We must perform the rotation independently for each Euler angle
% to accomodate having rotation centers that differ by Euler angle.
switch directionFlag
    case 'forward'
        
        % Torsion
        eyeCoord = eyeCoord - rotationCenters.tor;
        eyeCoord = (R.tor*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.tor;

        % Elevation
        eyeCoord = eyeCoord - rotationCenters.ele;
        eyeCoord = (R.ele*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.ele;

        % Azimuth
        eyeCoord = eyeCoord - rotationCenters.azi;
        eyeCoord = (R.azi*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.azi;                
                
    case 'inverse'
        
        % Azimuth
        eyeCoord = eyeCoord - rotationCenters.azi;
        eyeCoord = (R.azi'*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.azi;

        % Elevation
        eyeCoord = eyeCoord - rotationCenters.ele;
        eyeCoord = (R.ele'*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.ele;
        
        % Torsion
        eyeCoord = eyeCoord - rotationCenters.tor;
        eyeCoord = (R.tor'*eyeCoord')';
        eyeCoord = eyeCoord + rotationCenters.tor;
                        
end


end