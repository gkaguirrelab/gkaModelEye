

figure('NumberTitle', 'off', 'Name', 'Emmetropia')


% Create the eye and plot it
subplot(3,1,1)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{0}];
[navarroD, ~, path1, path2] = calcAccommodation(0,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','green');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','green','viewAngle',[0 90]);
xlim([-25 67]);
title('Focused at infinity (0D)')

% Now the accomodated eye
subplot(3,1,2)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{0}];
[navarroD, ~, path1, path2] = calcAccommodation(15,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','red');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','red','viewAngle',[0 90]);
xlim([-25 67]);
title('Focused at 67mm (15D)')

% Now the accomodated eye
subplot(3,1,3)
eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{15}];
[navarroD, ~, path1, path2] = calcAccommodation(15,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','green');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','green','viewAngle',[0 90]);
xlim([-25 67]);
title('Focused at 67mm (15D) using a spectacle lens')
