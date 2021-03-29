function plotHandles = addCameraIcon(sceneGeometry,aziEle)
% Adds a camera (with LED light sources) to the current image
%
% Syntax:
%  plotHandles = addCameraIcon(sceneGeometry)
%
% Description:
%   Renders a simple camera illustration in the currently active figure,
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
    plotOpticalSystem(sceneGeometry);
    addCameraIcon(sceneGeometry);
%}



if nargin==1
    aziEle = [0, 0];
end

%% Grab the nodal point of the camera
nodalPoint = convertWorldToEyeCoord(sceneGeometry.cameraPosition.translation);
sensorRatio = sceneGeometry.cameraIntrinsic.sensorResolution(2) ./ ...
    sceneGeometry.cameraIntrinsic.sensorResolution(1);

%% Define an empty variable to hold the plot objects
plotHandles = gobjects(0);

%% Plot the camera body
edges = [10, 40, 40*sensorRatio];
alpha = 0.5;
clr = [0.75 0.75 1];

% Figure out where the origin of the cuboid goes
trans = edges.*[0 0.5 0.5];
origin = nodalPoint - trans;

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


%% Plot the lens
S = quadric.scale(quadric.unitSphere,[2000 5 5]);
S = quadric.translate(S,nodalPoint);
boundingBox = [nodalPoint(1)-15 nodalPoint(1) nodalPoint(2)-5 nodalPoint(2)+5 nodalPoint(3)-5 nodalPoint(3)+5];
plotHandles(end+1) = quadric.plotGridSurface(S, boundingBox, [0.5 0.5 1], 0.5);

% Add some circles to cap the lens ends
boundingBox = [nodalPoint(1)-0.1 nodalPoint(1)+0.1 nodalPoint(2)-5 nodalPoint(2)+5 nodalPoint(3)-5 nodalPoint(3)+5];
plotHandles(end+1) = quadric.plotGridSurface(S, boundingBox, 'k', 0.5);
S = quadric.translate(S,[-15 0 0]);
boundingBox = [nodalPoint(1)-14.9 nodalPoint(1)-15.1 nodalPoint(2)-6 nodalPoint(2)+6 nodalPoint(3)-6 nodalPoint(3)+6];
plotHandles(end+1) = quadric.plotGridSurface(S, boundingBox, 'k', 0.5);

%% Add the glints
glintPoints = convertWorldToEyeCoord(sceneGeometry.cameraPosition.glintSourceRelative);
for gg = 1:size(glintPoints,1)
    coord = nodalPoint+glintPoints(gg,:);
    S = quadric.scale(quadric.unitSphere,[1 1 1]);
    S = quadric.translate(S,coord);
    boundingBox = [coord(1)-1 coord(1) coord(2)-1 coord(2)+1 coord(3)-1 coord(3)+1];
    plotHandles(end+1) = quadric.plotGridSurface(S, boundingBox, [1 0 0], 1);
end


%% Rotate the camera around the rotation center of the eye
rotate(plotHandles,[0 0 1],aziEle(1),sceneGeometry.eye.rotationCenters.azi);
rotate(plotHandles,[0 1 0],aziEle(2),sceneGeometry.eye.rotationCenters.ele);


end
