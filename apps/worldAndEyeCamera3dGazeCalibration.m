% worldAndEyeCamera3dGazeCalibration
%
% This script simulates the circumstance of a head-mounted light logger
% that includes a world camera and an eye camera. An application of such
% systems is to use the eye camera to estimate the location of gaze within
% images obtained from the world camera. We simulate here collection of 3D
% gaze calibration data for a set of targets

% First, create a scene geometry in which the coordinate system is
% organized around an eye camera that is aligned with the optical axis of
% the eye in primary position (i.e., straight ahead, The default option)
eye = modelEyeParameters();

% Calculate the eye rotation that would be needed to fixate upon the
% actual location of the infrared eye camera
eyeCameraTranslation = [-10,-15,25];

% Define the position in visual angle of the intended IR camera, w.r.t. the primary position of the
% eye
[fieldAngularPosition,targetDistance] = calcFieldAngularPosition(eye,eyeCameraTranslation);

% Obtain the eye pose needed to fixate this intended IR camera position
eyePose = calcFixationPose(eye,fieldAngularPosition,targetDistance);

% Now set up the sceneGeometry for the actual position of the IR camera. We
% will define the primary position of the eye as the negative of the eye
% pose calculated earlier.
primaryPosition = -eyePose(1:2);
cameraPosition.translation = eyeCameraTranslation;
cameraPosition.torsion = 0;
cameraPosition.glintSourceRelative = [1,0,1];
sceneGeometry = createSceneGeometry('cameraPosition',cameraPosition,'primaryPosition',primaryPosition);

% Define a set of gaze target that is at the 
gazeTargets = {...
    [0,0,1000]...
    };

% For each gaze target, find the fixation pose of the eye
