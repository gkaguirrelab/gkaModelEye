

figure('NumberTitle', 'off', 'Name', 'Emmetropia')

eyeVarargin = [{'sphericalAmetropia'},{0},{'spectacleLens'},{0}];

% Obtain the rays and the D parameter
[navarroD, ~, path1, path2] = calcAccommodation(0,eyeVarargin{:});

% Create the eye and plot it
subplot(2,1,1)
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);

% Add the rays
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','green');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','green','viewAngle',[0 90]);
xlim([-30 30]);

% Now the accomodated eye
subplot(2,1,2)
[navarroD, ~, path1, path2] = calcAccommodation(5,eyeVarargin{:});
sceneGeometry = createSceneGeometry('navarroD',navarroD,eyeVarargin{:});
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);

% Add the rays
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','red');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','red','viewAngle',[0 90]);
xlim([-30 30]);

