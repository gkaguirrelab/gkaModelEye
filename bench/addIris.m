function opticalSystemOut = addIris(opticalSystemIn, stopRadius, irisColor)
% Add an iris and stop to an eye optical system
%
% Syntax:
%   opticalSystemOut = addIris(opticalSystemIn, stopRadius)
%
% Description:
%   Adds the iris surface and aperture stop to an eye optical system
%   structure. The iris surface is not included by default in the optical
%   system as the radius of the stop is a dynamic parameter. The iris
%   surface is modeled as having an infinite index of refraction, so that
%   rays that intersect this surface go no further. The aperture stop is a
%   circulr, must-intersect surface.
%
%   If a iris is already present, then calling this routine will update the
%   surfaces with the passed stopRadius
%
% Inputs:
%   opticalSystemIn       - Struct. Must have the fields
%                             {opticalSystem, surfaceLabels, surfaceColors}
%                           See "assembleOpticalSystem.m" for details.
%   stopRadius            - Scalar. Radius of the stop in mm. Defaults to
%                           2 mm.
%   irisColor             - Char vector. Choose from 'brown', 'green', and
%                           'blue'.
%
% Outputs:
%   opticalSystemOut      - Struct.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    opticalSystemStruct = sceneGeometry.refraction.retinaToCamera;
    opticalSystemStruct = addIris(opticalSystemStruct);
    plotOpticalSystem('surfaceSet',opticalSystemStruct,'addLighting',true,'viewAngle',[0 90]);
%}


% Require that the opticalSystem is a full opticalSystem structure
if ~isstruct(opticalSystemIn)
    error('Please provide a full opticalSystem structure for this function')
end


%% Handle nargin
if nargin==1
    % Assume 2 mm
    stopRadius = 2;
    irisColor = 'brown';
end
if nargin==2
    irisColor = 'brown';
end


%% Set the color of the iris
% The 4th value is the alpha for the iris. It has a value of greater than
% 2, so that it will still be quite opaque once the optical system is
% rendered with an overal surface alpha of (e.g.) 0.1.
switch irisColor
    case 'brown'
        surfaceColor = [1.0000    0.7880    0.5252 2];
    case 'blue'
        surfaceColor = [0.1804    0.3255    1 2];
    case 'green'
        surfaceColor = [0.2392    1    0.1137 2];
    otherwise
        % Default to brown
        surfaceColor = [0.3882    0.3059    0.2039 2];
end


%% Place the opticalSystem in the eyeToCamera direction
if strcmp(calcSystemDirection(opticalSystemIn.opticalSystem),'cameraToEye')
    opticalSystemOut = reverseSystemDirection(opticalSystemIn);
    reverseSystemFlag = true;
else
    opticalSystemOut = opticalSystemIn;
    reverseSystemFlag = false;
end


%% Prep the optical system matrix
% Extract the opticalSystem matrix
opticalSystemMatrix = opticalSystemIn.opticalSystem;

% Strip the optical system of nan rows. We will add these back before
% returning
numRows = size(opticalSystemMatrix,1);
opticalSystemMatrix = opticalSystemMatrix(sum(isnan(opticalSystemMatrix),2)~=size(opticalSystemMatrix,2),:);


%% Set the target row
% Find the row of the optical system where we will place the iris.surface
% and iris.stop
irisPresent = any(strcmp(opticalSystemIn.surfaceLabels,'iris.surface'));
if irisPresent
    % We will update the existing surfaces
    targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'iris.surface'));
    replaceExistingFlag = true;
    surfaceColor = opticalSystemIn.surfaceColors{targetRow};
else
    targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'lens.front'));
    replaceExistingFlag = false;
end



%% Assemble the iris surface
eye.meta.eyeLaterality = 'Right';
iris = human.iris( eye );
irisCenterX = iris.center(1);

% The stop is an ellipsoid that is elongated, so that it has a flat
% surface as seen by rays traveling parallel to the optical axis. The stop
% has a trivial non-zero depth.
eRadii = [3 10 10];
S = quadric.scale(quadric.unitSphere,eRadii);

% Find the x axis position on the quadric where the radius in the y-z plane
% is equal to the iris radius, and the position where the radius is equal
% to the stop radius
F = quadric.vecToFunc(S);
myObj = @(x) abs(F(x,iris.radius,0));
xPerim = fminbnd(myObj,0,eRadii(1));
myObj = @(x) abs(F(x,stopRadius,0));
xStop = fminbnd(myObj,0,eRadii(1));

% Translate the quadric so that this xStop position corresponds to the iris
% center on the x axis
S = quadric.translate(S,[irisCenterX-xStop 0 0]);

% Assemble a line for the optical system
irisSurfaceLine = nan(1,19);

% Add the quadric itself
irisSurfaceLine(1:10) = quadric.matrixToVec(S);

% The system will be in eyeToCamera orientation, so rays traveling
% left-to-right encounter concave surface
irisSurfaceLine(11) = 1;

% Add the bounding box
bb_surface = [irisCenterX-xStop+xPerim irisCenterX -iris.radius iris.radius -iris.radius iris.radius];
irisSurfaceLine(12:17) = bb_surface;

% Now the stop
bb_stop = [irisCenterX irisCenterX+eRadii(1)-xStop -stopRadius stopRadius -stopRadius stopRadius];
irisStopLine = irisSurfaceLine;
irisStopLine(12:17) = bb_stop;

% This is the iris surface, so we don't want to intersect
irisSurfaceLine(18) = 0;

% This is the iris stop, so we do want to intersect
irisStopLine(18) = 1;

% Set the index of refraction of the surface to inf, causing the ray to
% stop if it hits the iris surface
irisSurfaceLine(19) = Inf;

% For the stop, use the same index of refraction as the surrounding medium,
% so that the aperture stop does not have any effect of refraction.
irisStopLine(19) = opticalSystemMatrix(targetRow,19);


%% Add the information to the optical system
surfaceLabels = opticalSystemIn.surfaceLabels;
surfaceColors = opticalSystemIn.surfaceColors;

if replaceExistingFlag
    opticalSystemMatrix(targetRow,:) = irisSurfaceLine;
    opticalSystemMatrix(targetRow+1,:) = irisStopLine;
else
    % Matrix
    opticalSystemMatrix = [ ...
        opticalSystemMatrix(1:targetRow,:); ...
        irisSurfaceLine; ...
        irisStopLine; ...
        opticalSystemMatrix(targetRow+1:end,:)];
    % Labels
    surfaceLabels = [surfaceLabels(1:targetRow); ...
        'iris.surface'; ...
        'iris.stop'; ...
        surfaceLabels(targetRow+1:end)];
    % Colors
    surfaceColors = [surfaceColors(1:targetRow); ...
        surfaceColor; ...
        [1 0 0 1]; ...
        surfaceColors(targetRow+1:end)];
end


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
if numRows > size(opticalSystemMatrix,1)
    opticalSystemMatrix = [opticalSystemMatrix; ...
        nan(numRows-size(opticalSystemMatrix,1),19)];
end


%% Handle labels, colors, and the output variable
% Handle the output variable
opticalSystemOut.opticalSystem = opticalSystemMatrix;
opticalSystemOut.surfaceLabels = surfaceLabels;
opticalSystemOut.surfaceColors = surfaceColors;


%% Reverse the system if needed
if reverseSystemFlag
    opticalSystemOut = reverseSystemDirection(opticalSystemOut);
end

end
