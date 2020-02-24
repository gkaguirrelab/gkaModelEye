function eyeCoord = rotateEyeCoord(eyeCoord, eyePose, rotationCenters)
% Apply an eye rotation to an eye coordinate
%
% Syntax:
%  eyeCoord = rotateEyeCoord(eyeCoord, eyePose, rotationCenters)
%
% Description
%   The eye coordinate space is defined along the optical axis of the eye
%   when the eye is aligned with the optical axis of a camera. This routine
%   takes an eye coordinate and updates the value by applying an eye
%   rotation.
%

% Generate the rotation matrix
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Apply the rotation, making use of the rotationCenters of the eye

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