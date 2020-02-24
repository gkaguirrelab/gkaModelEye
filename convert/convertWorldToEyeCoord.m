function eyeCoord = convertWorldToEyeCoord(eyePose, worldCoord, rotationCenters)
% World to eye coordinate conversion, accounting for eye rotation

% We assign the variable eyeCoordTarget the coordinates of the target after
% conversion to eye coordinates. Then, the point is counter-rotated by the
% eye pose, so that the eyeCoordTarget is in a position equivalent to if
% the eye had rotated.
cameraRot = -eyePose;
RotAzi = [cosd(cameraRot(1)) -sind(cameraRot(1)) 0; sind(cameraRot(1)) cosd(cameraRot(1)) 0; 0 0 1];
RotEle = [cosd(-cameraRot(2)) 0 sind(-cameraRot(2)); 0 1 0; -sind(-cameraRot(2)) 0 cosd(-cameraRot(2))];
RotTor = [1 0 0; 0 cosd(cameraRot(3)) -sind(cameraRot(3)); 0 sind(cameraRot(3)) cosd(cameraRot(3))];

% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. This is now the eyeTarget or ET
eyeCoord = worldCoord([3 1 2])';

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