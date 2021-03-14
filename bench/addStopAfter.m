function opticalSystemOut = addStopAfter(opticalSystemIn, stopRadius, stopCenter, targetRow, replaceExistingFlag, surfaceLabel, surfaceColor)
% Add a stop after the specified surface position
%
% Syntax:
%   opticalSystemOut = addStopAfter(opticalSystemIn, stopRadius, stopCenter, targetRow, replaceExistingFlag, surfaceLabel, surfaceColor)
%
% Description:
%   Adds a stop to an optical system centered on the optical axis. The stop
%   is a must-intersect, circular surface of the specified radius and
%   position. The optical system can be provided as a full structure
%   (including labels and plot colors) or as a matrix. If an opticalSystem
%   structure is provide with minimal other inputs,the routine will assume
%   that the stop is the iris stop of the eye, and will update the radius
%   of this stop if it already exists.
%
% Inputs:
%   opticalSystemIn       - Struct or matrix. If struct, must have the
%                           fields {opticalSystem, surfaceLabels,
%                           surfaceColors}. The matrix form is mx19. See
%                           "assembleOpticalSystem.m" for details. 
%   stopRadius            - Scalar. Radius of the stop in mm. Defaults to
%                           2 mm.
%   stopCenter            - Scalar. The position of the stop on the optical 
%                           axis in mm. Defaults to -3.9 mm.
%   targetRow             - Scalar. The row of the optical system after
%                           which the stop is to be inserted, or the row
%                           that contains the surface to be replaced.
%   replaceExistingFlag   - Logical. Controls if the stop is to be added to
%                           the opticalSystem or replace an existing row.
%   surfaceLabel          - Char vector or string. The label for the stop.
%   surfaceColor          - 3x1 vector that specifies the color to be used
%                           to display the stop in the plotOpticalSystem
%                           function.
%
% Outputs:
%   opticalSystemOut      - Struct or matrix. Will match the type of
%                           opticalSystemIn
%
% Examples:
%{
    sceneGeometry = createSceneGeometry('spectralDomain','vis');
    opticalSystemStruct = sceneGeometry.refraction.retinaToCamera;
    opticalSystemStruct = addStopAfter(opticalSystemStruct);
    plotOpticalSystem(opticalSystemStruct,'viewAngle',[90 0]);
%}


%% Handle nargin
switch nargin
    case 1
        % Assume we are adding a 2mm iris aperture stop. Requires that the
        % opticalSystemIn be a full opticalSystem structure.
        if isstruct(opticalSystemIn)
            % Assume 2 mm
            stopRadius = 2;
            % Assign a stop center location that is the model location for
            % the iris
            stopCenter = [-3.9 0 0];
            % Check if there is an existing iris.stop surface
            irisStopPresent = any(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
            if irisStopPresent
                % We will update the existing iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
                replaceExistingFlag = true;
                surfaceLabel = 'iris.stop';
                surfaceColor = opticalSystemIn.surfaceColors{targetRow};
            else
                % Find the row of the optical system that contains the
                % front of the lens and we will add an iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'lens.front'));
                replaceExistingFlag = false;
                surfaceLabel = 'iris.stop';
                surfaceColor = [0 0 0];
            end
        else
            error('gkaModelEye:addStopAfter','addStopAfter requires at least an opticalSystem structure as input');
        end
    case 2
        % Assume we are adding an iris aperture stop of specified radius.
        % Requires that the opticalSystemIn be a full opticalSystem
        % structure.
        if isstruct(opticalSystemIn)
            % Assign a stop center location that is the model location for
            % the iris
            stopCenter = [-3.9 0 0];
            % Check if there is an existing iris.stop surface
            irisStopPresent = any(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
            if irisStopPresent
                % We will update the existing iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
                replaceExistingFlag = true;
                surfaceLabel = 'iris.stop';
                surfaceColor = opticalSystemIn.surfaceColors{targetRow};
            else
                % Find the row of the optical system that contains the
                % front of the lens and we will add an iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'lens.front'));
                replaceExistingFlag = false;
                surfaceLabel = 'iris.stop';
                surfaceColor = [0 0 0];
            end
        else
            error('gkaModelEye:addStopAfter','addStopAfter requires at least an opticalSystem structure as input');
        end
    case 3
        % Assume we are adding an iris aperture stop of specified radius
        % and center location.  Requires that the opticalSystemIn be a full
        % opticalSystem structure.
        if isstruct(opticalSystemIn)
            % Check if there is an existing iris.stop surface
            irisStopPresent = any(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
            if irisStopPresent
                % We will update the existing iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'iris.stop'));
                replaceExistingFlag = true;
                surfaceLabel = 'iris.stop';
                surfaceColor = opticalSystemIn.surfaceColors{targetRow};
            else
                % Find the row of the optical system that contains the
                % front of the lens and we will add an iris.stop
                targetRow = find(strcmp(opticalSystemIn.surfaceLabels,'lens.front'));
                replaceExistingFlag = false;
                surfaceLabel = 'iris.stop';
                surfaceColor = [0 0 0];
            end
        else
            error('gkaModelEye:addStopAfter','addStopAfter requires at least an opticalSystem structure as input');
        end
    case 4
        % We have sufficient inputs that we can accept an opticalSystem
        % matrix or an opticalSystem structure.
        replaceExistingFlag = false;
        surfaceLabel = 'stop';
        surfaceColor = [0 0 0];
    case 5
        surfaceLabel = 'stop';
        surfaceColor = [0 0 0];
    case 6
        surfaceLabel = 'stop';
    case 7
        % all is well
    otherwise
        error('gkaModelEye:addStopAfter','Invalid number of arguments.');
end

% Check if the stopCenter is a scalar, in which case make it a coordinate
if isscalar(stopCenter)
    stopCenter = [stopCenter 0 0];
end


%% Prepeare the optical system

% Extract the opticalSystem matrix
opticalSystemOut = opticalSystemIn;
if isstruct(opticalSystemIn)
    opticalSystemMatrix = opticalSystemIn.opticalSystem;
else
    opticalSystemMatrix = opticalSystemIn;
end

% Determine system direction
systemDirection = calcSystemDirection(opticalSystemMatrix);

% Strip the optical system of nan rows. We will add these back before
% returning
numRows = size(opticalSystemMatrix,1);
opticalSystemMatrix = opticalSystemMatrix(sum(isnan(opticalSystemMatrix),2)~=size(opticalSystemMatrix,2),:);


%% Assemble the stop surface

% The stop is an ellipsoid that is very elongated, so that it has a flat
% surface as seen by rays traveling parallel to the optical axis. The stop
% has a trivial non-zero depth.
ellipseDepth = 0.01;
S = quadric.scale(quadric.unitSphere,[ellipseDepth,stopRadius*5,stopRadius*5]);
t = stopCenter; t(1) = t(1)-ellipseDepth;
S = quadric.translate(S,t);
stopFront = stopCenter(1);

% Find the x-axis position at which the height of the ellipsoid is equal to
% the desired stop radius
F = quadric.vecToFunc(S);
myObj = @(x) radiusAtX(F,stopFront-ellipseDepth+x)-stopRadius;
x = fzero(myObj,ellipseDepth-(1e-4));

% Assemble a line for the optical system
stopLine = nan(1,19);

% Add the quadric itself
stopLine(1:10) = quadric.matrixToVec(S);

% The side of the stop that rays encounter varies for the system direction
switch systemDirection
    case 'eyeToCamera'
        % Rays traveling left-to-right encounter concave surface
        stopLine(11) = 1; 
    case 'cameraToEye'
        % Rays traveling left-to-right encounter convex surface
        stopLine(11) = -1; 
end

% Add the bounding box
bb = [stopFront-x stopFront  -stopRadius stopRadius -stopRadius stopRadius];
stopLine(12:17) = bb;

% The ray must intersect the stop
stopLine(18) = 1;

% Use the same index of refraction as the surrounding medium, so that the
% aperture stop does not have any effect of refraction.
stopLine(19) = opticalSystemMatrix(targetRow,19);

% Add or update the stop in the optical system matrix
if replaceExistingFlag
    % Replace the targetRow
    opticalSystemMatrix(targetRow,:) = stopLine;
else
    % Add this line to the opticalSystemMatrix after the targetRow
    opticalSystemMatrix = [ ...
        opticalSystemMatrix(1:targetRow,:); ...
        stopLine; ...
        opticalSystemMatrix(targetRow+1:end,:)];

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
if replaceExistingFlag
    
    % Handle the output variable
    if isstruct(opticalSystemIn)
        opticalSystemOut.opticalSystem = opticalSystemMatrix;
    else
        opticalSystemOut = opticalSystemMatrix;
    end
else
    
    % Handle the output variable
    if isstruct(opticalSystemIn)
        opticalSystemOut.opticalSystem = opticalSystemMatrix;

        % If we have an opticalSystem struct, add a label and plot colors        
        surfaceLabels = opticalSystemIn.surfaceLabels;
        surfaceLabels = [surfaceLabels(1:targetRow); ...
            surfaceLabel; ...
            surfaceLabels(targetRow+1:end)];
        surfaceColors = opticalSystemIn.surfaceColors;
        surfaceColors = [surfaceColors(1:targetRow); ...
            surfaceColor; ...
            surfaceColors(targetRow+1:end)];
        opticalSystemOut.surfaceLabels = surfaceLabels;
        opticalSystemOut.surfaceColors = surfaceColors;
    else
        opticalSystemOut = opticalSystemMatrix;
    end
end


end


%% LOCAL FUNCTION

function z = radiusAtX(F,x)
    myObj = @(z) F(x,0,z);
    z = abs(fzero(myObj,0));
end


