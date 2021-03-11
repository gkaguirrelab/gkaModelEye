%% generateImages
% This script creates and saves the images that are displayed in the README
% file for the repo.

% The location of this function, which is where we will save the images
functionDirPath = fileparts(mfilename('fullpath'));



%% Create the model to plot

% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another. Set the key-value to calculate the
% location of the fovea
sceneGeometry=createSceneGeometry();

% Calculate the line-of-sight axis of the eye
rayDestination = sceneGeometry.eye.landmarks.fovea.coords;
rayPathLoS = calcSightRayToRetina(sceneGeometry.eye,rayDestination);

% Add an iris and aperture stop for the optical system render
sceneGeometry.refraction.retinaToCamera = addIris(sceneGeometry.refraction.retinaToCamera, 2, 'green');


%% opticalSystem3D

% Set up a figure
figHandle = figure('Visible','off');
set(gcf,'PaperOrientation','portrait');
set(figHandle, 'Units','inches')
height = 4;
width = 4;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');

% Coordinates of the fovea
p = sceneGeometry.eye.landmarks.fovea.coords';

% Assemble the ray
R = quadric.normalizeRay( ...
    quadric.anglesToRay(p, ...
    sceneGeometry.eye.landmarks.fovea.degField(1), ...
    sceneGeometry.eye.landmarks.fovea.degField(2)));

% Add this ray to the optical system plot
plotOpticalSystem(sceneGeometry,'surfaceSet','retinaToCamera','newFigure',false,'rayPath',rayPathLoS);
xlim([-25 5]);

% Save this image
filename = fullfile(functionDirPath,'opticalSystem3D.png');
print(figHandle,filename,'-dpng','-r300');

% Close the figure
close(figHandle)


%% renderEyePose

% Define the eyePose and a background image
eyePose = [20 5 0 3];
backgroundImage = zeros(480,640)+0.5;

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y' 'Qr'};

% Render
figHandle = renderEyePose(eyePose, sceneGeometry,...
    'newFigure',true,'visible',false,...
    'modelEyeLabelNames',modelEyeLabelNames,...
    'modelEyePlotColors',modelEyePlotColors,...
    'modelEyeSymbolSizeScaler',0.75,'backgroundImage',backgroundImage);

% Set up a figure
set(figHandle,'PaperOrientation','portrait');
set(figHandle, 'Units','inches')
height = 4;
width = 4;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color',[0.5 0.5 0.5],...
    'InvertHardCopy', 'off');

% Save this image
filename = fullfile(functionDirPath,'renderEyePose.png');
print(figHandle,filename,'-dpng','-r300');

% Close the figure
close(figHandle)


%% modelEyeSchematic

% Set up a figure
figHandle = figure('Visible','off');
set(gcf,'PaperOrientation','portrait');
set(figHandle, 'Units','inches')
height = 4;
width = 4;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');

% Plot the schematic eye in red
plotModelEyeSchematic(sceneGeometry.eye,'view','horizontal',...
    'newFigure',false,'plotColor','r', ...
    'rayPath',rayPathLoS);

% Remove the fovea
sceneGeometry.eye.landmarks=rmfield(sceneGeometry.eye.landmarks,'fovea');

% Now over-plot in black, without the fovea and line of sight
plotModelEyeSchematic(sceneGeometry.eye,'view','horizontal','newFigure',false,'plotColor','k')

% Clean up the plot limits
xlim([-25 5])
ylim([-15 15])
axis square

% Save this image
filename = fullfile(functionDirPath,'modelEyeSchematic.png');
print(figHandle,filename,'-dpng','-r300');

% Close the figure
close(figHandle)
