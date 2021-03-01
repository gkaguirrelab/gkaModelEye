%% DEMO_accommodation
% Demonstrate adjustment of accommodation of the model eye
%
% Description:
%   The accommodation of the eye is controled by passing values for the
%   navarroD key to the modelEyeParameters.m function. To find the navarroD
%   value that corresponds to a particular accommodative state of the eye,
%   use the function calcAccommodation.
%


% Create a figure
figure('NumberTitle', 'off', 'Name', 'Varying accommodative state')

% Create the eye and plot it
subplot(3,1,1)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{0}];
[navarroD, ~, path1, path2] = calcAccommodation(0,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','green');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','green','viewAngle',[0 90]);
ylim([-25 25]);
xlim([-25 67]);
title('Focused at infinity (0D)')
drawnow

% Now the accomodated eye
subplot(3,1,2)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{0}];
[navarroD, ~, path1, path2] = calcAccommodation(15,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','red');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','red','viewAngle',[0 90]);
ylim([-25 25]);
xlim([-25 67]);
title('Focused at 67mm (15D)')
drawnow

% Now the eye looking through a +15 magnifying glass
subplot(3,1,3)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{15}];
[navarroD, ~, path1, path2] = calcAccommodation(15,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','red');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','red','viewAngle',[0 90]);
ylim([-25 25]);
xlim([-25 67]);
title('Focused at 67mm (15D) using a magnifying glass (spectacle lens)')
