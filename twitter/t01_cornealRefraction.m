%% t01_cornealRefraction
% Demonstrate the effect of refraction by the cornea on the appearance of
% the entrance pupil of the eye


% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% Define an eyePose with the azimuth, elevation, torsion rotation of the 
% eye in degrees, and the radius of the aperture stop of the iris in mm
eyePose = [-30 -5 0 3];

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% Render the eye with corneal refraction
renderEyePose(eyePose, sceneGeometry, ...
    'modelEyeLabelNames',modelEyeLabelNames,...
    'modelEyePlotColors',modelEyePlotColors);

% Remove the effect of corneal refraction. To do so, we will hack the
% index of refraction of the aqueous humor, cornea, and tear film to all
% be one.
sceneGeometry.refraction.stopToMedium.opticalSystem(1:3,19) = 1;

% Indicate that we only want to render the "pupil", which is now the same
% as the aperture stop, as we have removed the refractive effects of the
% cornea.
modelEyeLabelNames = {'pupilEllipse'};
modelEyePlotColors = {'-r'};

% Render the eye again
renderEyePose(eyePose, sceneGeometry, ...
    'newFigure', false,...
    'modelEyeLabelNames',modelEyeLabelNames,...
    'modelEyePlotColors',modelEyePlotColors);

% Add a title.
title('Stop in red, pupil in green')
