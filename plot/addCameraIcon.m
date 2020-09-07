function plotHandles = addCameraIcon(sceneGeometry)
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
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    addCameraIcon(sceneGeometry);
%}


%% Grab the nodal point of the camera
nodalPoint = convertWorldToEyeCoord(sceneGeometry.cameraPosition.translation);


%% Define an empty variable to hold the plot objects
plotHandles = gobjects(0);

%% Plot the camera body
edges = [10, 40, 20];
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
plotHandles(end+1) = quadric.plotSurface(S, boundingBox, [0.5 0.5 1], 0.5);

% Add some circles to cap the lens ends
boundingBox = [nodalPoint(1)-1 nodalPoint(1)+1 nodalPoint(2)-5 nodalPoint(2)+5 nodalPoint(3)-5 nodalPoint(3)+5];
plotHandles(end+1) = quadric.plotSurface(S, boundingBox, 'none', 0, 'none', 'black');
S = quadric.translate(S,[-15 0 0]);
boundingBox = [nodalPoint(1)-16 nodalPoint(1)+14 nodalPoint(2)-6 nodalPoint(2)+6 nodalPoint(3)-6 nodalPoint(3)+6];
plotHandles(end+1) = quadric.plotSurface(S, boundingBox, 'none', 0, 'none', 'black');


%% Add the glints
glintPoints = convertWorldToEyeCoord(sceneGeometry.cameraPosition.glintSourceRelative);
for gg = 1:size(glintPoints,1)
    coord = nodalPoint+glintPoints(gg,:);
    S = quadric.scale(quadric.unitSphere,[1 1 1]);
    S = quadric.translate(S,coord);
    boundingBox = [coord(1)-1 coord(1) coord(2)-1 coord(2)+1 coord(3)-1 coord(3)+1];
    plotHandles(end+1) = quadric.plotSurface(S, boundingBox, [1 0 0], 1);
end

end
