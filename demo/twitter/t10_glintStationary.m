%% t10_glintStationary
% Illustrate the ray path from light to camera to produce a glint
%{
10: The "glint" is the reflection of a light source from the front corneal surface (the 1st Purkinje image). Shown is the ray path from an active IR source that reflects off the eye and returns to the stop of a camera.
%}
% For details see:
%   project/stages/addGlint.m
%

% Create the default scene geometry
sceneGeometry = createSceneGeometry();

% Define an eye pose
eyePose = [0 0 0 2];

% The position of the glint source in world coordinates
glintSourceWorld = sceneGeometry.cameraPosition.translation + ...
    sceneGeometry.cameraPosition.glintSourceRelative;

% Assemble the args for the glint ray trace
args = {sceneGeometry.cameraPosition.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.refraction.cameraToMedium.opticalSystem, ...
    sceneGeometry.refraction.glint.opticalSystem, ...
    sceneGeometry.refraction.mediumToCamera.opticalSystem};

% Perform the computation using the passed function handle
[~, initialRay] = ...
    findGlintRay(glintSourceWorld, eyePose, args{:});

% Ray trace from the light source, to the tear film, and back to the camera
[outputRay, rayPath] = rayTraceQuadrics(initialRay', sceneGeometry.refraction.glint.opticalSystem);

% Create a figure
figure
set(gcf,'color','w');

% Render the optical system and ray path
plotOpticalSystem(sceneGeometry, ...
    'outputRay',outputRay,'rayPath',rayPath,...
    'outputRayScale',sceneGeometry.cameraPosition.translation(3),...
    'newFigure',false);

% Add the camera
addCameraIcon(sceneGeometry);

% Provide a title
title('Ray path from light source to camera')
