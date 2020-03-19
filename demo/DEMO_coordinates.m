%% Demo coordinate axes
% Plot and label the axes used in the code and projections
% 
% Description
%   The representation of the eye and the scene in the model code makes use
%   of two different coordinate frames. This script creates a labeled 3D
%   plot that illustrates these frames and their relationship.
%


% Housekeeping
clear

sceneGeometry = createSceneGeometry();

figure


subplot(1,3,1)

% Add the eye
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'newFigure',false,'surfaceAlpha',0.4);

% Add axis lines for the eyeWorld space
plot3([-25 25],[0 0],[0 0],'-k')
plot3([0 0],[-25 25],[0 0],'-k')
plot3([0 0],[0 0],[-25 25],'-k')

% Add end-axis labels
text(-27,0,0,'-p1','HorizontalAlignment','center'); text(+27,0,0,'+p1','HorizontalAlignment','center');
text(0,-27,0,'-p2','HorizontalAlignment','center'); text(0,+27,0,'+p2','HorizontalAlignment','center');
text(0,0,-27,'-p3','HorizontalAlignment','center'); text(0,0,+27,'+p3','HorizontalAlignment','center');

% Indicate floor
f = fill3([-35 35 -35 35],[35 35 -35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)
f = fill3([-35 -35 35 35],[35 -35 35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)

% Handle title and axis
title('eye coordinates')
axis off



subplot(1,3,2)

% Add the eye
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'newFigure',false,'surfaceAlpha',0.4);

% Add axis lines for the eyeWorld space
plot3([-25 25],[0 0],[0 0],'-r')
plot3([0 0],[-25 25],[0 0],'-r')
plot3([0 0],[0 0],[-25 25],'-r')

% Add end-axis labels
text(-27,0,0,'-z','HorizontalAlignment','center'); text(+27,0,0,'+z','HorizontalAlignment','center');
text(0,-27,0,'-x','HorizontalAlignment','center'); text(0,+27,0,'+x','HorizontalAlignment','center');
text(0,0,-27,'-y','HorizontalAlignment','center'); text(0,0,+27,'+y','HorizontalAlignment','center');

% Indicate floor
f = fill3([-35 35 -35 35],[35 35 -35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)
f = fill3([-35 -35 35 35],[35 -35 35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)

% Handle title and axis
title('world coordinates')
axis off


subplot(1,3,3)

% Add the eye
imageSizeX = sceneGeometry.cameraIntrinsic.sensorResolution(1);
imageSizeY = sceneGeometry.cameraIntrinsic.sensorResolution(2);
backgroundImage = zeros(imageSizeX,imageSizeX)+0.75;
imshow(backgroundImage,[], 'Border', 'tight');
hold on

% Translate the camera up a bit so that the rendered eye is drawn in the
% middle of the frame
sceneGeometry.cameraPosition.translation(2)=4;

backgroundImage = zeros(imageSizeY,imageSizeX)+0.75;
eyePose = [0 0 0 3];
renderEyePose(eyePose, sceneGeometry,...
    'newFigure',false,...
    'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'},...
	'modelEyePlotColors', {'.w' 'ob' '-g' '.y'},...
	'modelEyeSymbolSizeScaler',0.5);

hold on

% Add some coordinate plots
text(50,40,'[0 0]','HorizontalAlignment','center');
text(400,40,'+x','HorizontalAlignment','left');
text(50,400,'+y','HorizontalAlignment','center');
plot([150 375],[40 40],'-b')
plot([50 50],[100 375],'-b')

% Handle title and axis
ylim([-25 680.5]);
title('image coordinates')
