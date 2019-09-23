

figure('NumberTitle', 'off', 'Name', 'Emmetropia')

% Create the eye and plot it
sceneGeometry = createSceneGeometry('sphericalAmetropia',0,'accommodationDiopters',0);
opticalSystem = sceneGeometry.refraction.cameraToRetina.opticalSystem;
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);

% Create parallel rays in the valid direction
R1 = quadric.normalizeRay([66.6667,-1;-0.3,0;0,0]);
R2 = quadric.normalizeRay([66.6667,-1;0.3,0;0,0]);

% Trace the rays
[~,path1] = rayTraceQuadrics(R1, opticalSystem);
[~,path2] = rayTraceQuadrics(R2, opticalSystem);

% Add the rays
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','green');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','green','viewAngle',[0 90]);

% Create the eye; no need to replot it as there is minimal visible change
% in the lens with the 
sceneGeometry = createSceneGeometry('sphericalAmetropia',0,'accommodationDiopters',15);
opticalSystem = sceneGeometry.refraction.cameraToRetina.opticalSystem;
%plotOpticalSystem('surfaceSet',sceneGeometry.refraction.cameraToRetina,'addLighting',true);

% Create diverging rays from the focal point
R1 = quadric.normalizeRay([66.6667,-1;0,-0.025;0,0]);
R2 = quadric.normalizeRay([66.6667,-1;0,+0.025;0,0]);

% Trace the rays
[~,path1] = rayTraceQuadrics(R1, opticalSystem);
[~,path2] = rayTraceQuadrics(R2, opticalSystem);

% Add the rays
plotOpticalSystem('newFigure',false,'rayPath',path1,'rayColor','red');
plotOpticalSystem('newFigure',false,'rayPath',path2,'rayColor','red','viewAngle',[0 90]);
