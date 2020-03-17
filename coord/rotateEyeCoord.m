function [eyeCoord, R] = rotateEyeCoord(eyeCoord, eyePose, rotationCenters, directionFlag)
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
end

%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to rightward, upward and clockwise (as seen
% from the perspective of the subject) eye movements
R.azi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
R.ele = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
R.tor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

switch directionFlag
    case 'forward'
        % This order (tor-ele-azi) corresponds to a head-fixed, extrinsic,
        % rotation matrix. The reverse order (azi-ele-tor) would be an
        % eye-fixed, intrinsic rotation matrix and would corresponds to the
        % "Fick coordinate" scheme.
        rotOrder = {'tor','ele','azi'};
        for rr = 1:3
            % We shift the points to each rotation center, rotate, shift
            % back, and repeat. We must perform the rotation independently
            % for each Euler angle to accomodate having rotation centers
            % that differ by Euler angle.
            eyeCoord = eyeCoord - rotationCenters.(rotOrder{rr});
            eyeCoord = (R.(rotOrder{rr})*eyeCoord')';
            eyeCoord = eyeCoord + rotationCenters.(rotOrder{rr});
            
        end
    case 'inverse'
        % Conduct the inverse rotation
        rotOrder = {'azi','ele','tor'};
        for rr = 1:3
            eyeCoord = eyeCoord - rotationCenters.(rotOrder{rr});
            eyeCoord = (R.(rotOrder{rr})'*eyeCoord')';
            eyeCoord = eyeCoord + rotationCenters.(rotOrder{rr});
        end
end


end