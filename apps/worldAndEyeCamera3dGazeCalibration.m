% worldAndEyeCamera3dGazeCalibration
%
% This script simulates the circumstance of a head-mounted light logger
% that includes a world camera and an eye camera. An application of such
% systems is to use the eye camera to estimate the location of gaze within
% images obtained from the world camera. We simulate here collection of 3D
% gaze calibration data for a set of targets. This simulation demonstrates
% that parallax errors are present in trying to relate eye movements to
% image collected from a world camera.

% House keeping
close all
clear

% Create the sceneGeometry
cameraTranslation = [-20;-10;100];
cameraGlintSourceRelative = [1;0;1];
eyeCameraGeometry = createSceneGeometry(...
    'cameraTranslation',cameraTranslation,...
    'cameraGlintSourceRelative',cameraGlintSourceRelative);

% Define a set of gaze targets
worldPoints = [...
    0,0,1000;...
    50,0,1000;...
    0,50,1000;...
    -50,0,1000;...
    0,-50,1000;...
    50,50,1000;...
    -50,-50,1000;...
    -50,50,1000;...
    50,-50,1000;...
        0,0,1500;...
    50,0,1500;...
    0,50,1500;...
    -50,0,1500;...
    0,-50,1500;...
    50,50,1500;...
    -50,-50,1500;...
    -50,50,1500;...
    50,-50,1500];

% For each gaze target, find the fixation pose of the eye
for gg = 1:size(worldPoints,1)
    targetWorldCoordinate = worldPoints(gg,:)';
    [fieldAngularPosition,targetDistance] = ...
        calcFieldAngularPosition(eyeCameraGeometry.eye,targetWorldCoordinate);
    eyePoses(gg,:) = calcFixationPose(eyeCameraGeometry.eye,fieldAngularPosition,targetDistance);
end

% Create an image of the targets from the perspective of the world camera
worldCameraTranslation = [-10,10,0];
worldCameraRotation = [0;0;0];
worldCameraGeometry = createSceneGeometry(...
    'cameraTranslation',cameraTranslation,...
    'cameraRotation',worldCameraRotation);
imagePoints = projectToImagePlane(worldPoints,worldCameraGeometry);
imagePoints = applyRadialLensDistortion(imagePoints,worldCameraGeometry);

% Grab the image size
imageSizeX = worldCameraGeometry.cameraIntrinsic.sensorResolution(1);
imageSizeY = worldCameraGeometry.cameraIntrinsic.sensorResolution(2);

% A blank frame to initialize the figure
backgroundImage = zeros(imageSizeY,imageSizeX)+0.5;

% Show the appearance of the targets in the world camera image
figure
imshow(backgroundImage,[], 'Border', 'tight');
axis off
axis equal
xlim([0 imageSizeX]);
ylim([0 imageSizeY]);
hold on
for idx = 1:size(imagePoints,1)
    markerSize = 5e4/worldPoints(idx,3);
    scatter(imagePoints(idx,1), imagePoints(idx,2), markerSize, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none');
    hold on
end
title('Targets in world camera image')
xlabel('horizontal [pixels]');
ylabel('vertical [pixels]');

% Plot the eye rotations
figure
for idx=1:size(eyePoses,1)
    markerSize = 5e4/worldPoints(idx,3);
    scatter(eyePoses(idx,1), eyePoses(idx,2), markerSize, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none');
    hold on
end
title('Eye rotations to targets')
xlabel('azimuth [deg]');
ylabel('elevation [deg]');
