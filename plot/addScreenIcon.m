function plotHandles = addScreenIcon(sceneGeometry,aziEle, gazeTargets)
% Adds a screen (with targets) to the current image
%
% Syntax:
%  plotHandles = addCameraIcon(sceneGeometry)
%
% Description:
%   Renders a simple screen illustration in the currently active figure,
%   with the position of the camera and LEDs taken from the passed
%   sceneGeometry.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Outputs:
%   plotHandles           - Array of handles for each of the surfaces
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    addCameraIcon(sceneGeometry);
%}

if nargin==1
    aziEle = [0, 0];
    gazeTargets = [];
end

if nargin==2
    gazeTargets = [];
end


%% Grab the nodal point of the camera
screenCenter = [sceneGeometry.screenPosition.screenDistance, 0, 0];


%% Define an empty variable to hold the plot objects
plotHandles = gobjects(0);


%% Plot the screen
edges = [1, sceneGeometry.screenPosition.dimensions(1), sceneGeometry.screenPosition.dimensions(2)];
alpha = 0.25;
clr = [0.5 0.5 0.5];

% Figure out where the origin of the cuboid goes
trans = edges.*[0 0.5 0.5];
origin = screenCenter - trans;

XYZ = { ...
    [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
    [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
    [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
    [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
    [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
    [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
    };
XYZ = mat2cell(...
    cellfun( @(x,y,z) x*y+z , ...
    XYZ , ...
    repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
    repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
    'UniformOutput',false), ...
    6,[1 1 1]);
plotHandles(end+1:end+6) = cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
    repmat({clr},6,1),...
    repmat({'FaceAlpha'},6,1),...
    repmat({alpha},6,1)...
    );

% If we have been given gazeTargets convert these from degrees to locations
% on the screen and plot them
if ~isempty(gazeTargets)
    
    % Visual angle turns out to be a surprisingly nuanced topic. For now,
    % we will calculate visual angle w.r.t. the optical center of the eye.
    opticalCenter = calcOpticalCenter(sceneGeometry.eye);
     
    % For each of the points, place a target at the center of the screen,
    % and then rotate it around the opticalCenter by the specied gaze
    % angle. I probably should introduce some pseudo torsion here to
    % correct for the Fick coordinate approach I am taking, but since this
    % is just for display purposes I am willing to let it go.
    for gg = 1:size(gazeTargets,2)
        plotHandles(end+1) = plot3(screenCenter(1),screenCenter(2),screenCenter(3),'.r','MarkerSize',30);
        rotate(plotHandles(end),[0 0 1],gazeTargets(1,gg),opticalCenter);
        rotate(plotHandles(end),[0 1 0],gazeTargets(2,gg),opticalCenter);
    end
end


%% Rotate the screen around the rotation center of the eye
rotate(plotHandles,[0 0 1],aziEle(1),sceneGeometry.eye.rotationCenters.azi);
rotate(plotHandles,[0 1 0],aziEle(2),sceneGeometry.eye.rotationCenters.ele);


end
