function [figHandle, plotHandles] = plotOpticalSystem(opticalSystem,options)
% Create a 3D rendered plot of the specified optical system
%
% Syntax:
%  figHandle = plotOpticalSystem(varargin)
%
% Description:
%   Create a 3D rendered plot of the specified optical system
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'newFigure'            - Logical. Determines if we create a new figure,
%                           or plot in the currently active figure.
%  'visible'              - Logical. If we are crearing a new figure,
%                           dertermines if that figure is made visible or
%                           not.
%  'surfaceSet'           - Struct or numeric. This is the specification of
%                           the optical system to be shown. If the calling
%                           function the optical system within a
%                           sceneGeometry struct, then this key value
%                           should be one of the fields within
%                               sceneGeometry.refraction
%                           with typical options being 'mediumToRetina' or
%                           'stopToMedium'. In this circumstance, the
%                           passed structure contains the fields
%                               opticalSystem, surfaceLabels,surfaceColors
%                           with the last two fields defining the plotted
%                           appearance of the surfaces. See
%                               model/createSceneGeometry.m
%                           for more information on the sceneGeometry
%                           structure. The routine also accepts just an
%                           opticalSystem matrix, in which case the system
%                           will be plotted in gray. See
%                               quadric/rayTraceQuadrics.m
%                           for a description of the opticalSystem matrix.
%  'surfaceAlpha'         - Scalar. This value scales the overall
%                           transparency of the rendering.
%  'retinaGeodetics'      - Logical. If set to true, and the surfaceSet is
%                           a structure that contains a "retina"
%                           surfaceLabel, then geodetic lines will be
%                           added to the retinal surface.
%  'rayPath'              - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(:,1)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan. This rayPath could be obtained
%                           using the function:
%                               quadric/rayTraceQuadrics.m
%  'rayColor'             - Char vector or 1x3 vector that specifies the
%                           line color for the ray.
%  'outputRay'            - 3x2 matrix that specifies the ray as a unit
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%  'outputRayColor'       - Char vector or 1x3 vector that specifies the
%                           line color for the outputRay.
%  'addLighting'          - Logical. Adds gouraud lighting to the plot.
%  'viewAngle'            - 1x2 vector. The view angle for the plot
%
% Outputs:
%   figHandle             - Handle to a created figure. Empty if a new
%                           figure was not requested.
%   plotHandles           - Array of handles for each of the surfaces
%
% Examples:
%{
    % Plot a selected surface set
    sceneGeometry = createSceneGeometry();
    plotOpticalSystem(sceneGeometry,'surfaceSet','retinaToCamera');
%}
%{
    %% Rays from the retina through the eye and a spectacle lens
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',-2);
    % Plot the optical system
    plotOpticalSystem(sceneGeometry,'surfaceSet','retinaToCamera');
    % Obtain the line-of-sight ray
    rayPath = calcSightRayToRetina(sceneGeometry.eye);
    % Add this ray to the optical system plot
    plotOpticalSystem('newFigure',false,'rayPath',rayPath);
    xlim([-30 30]);
%}


arguments
    opticalSystem {mustBeOpticalSystemCapable} = []
    options.newFigure (1,1) logical = true
    options.visible (1,1) logical = true
    options.surfaceSet char = 'mediumToRetina'
    options.surfaceAlpha (1,1) {mustBeNumeric} = 0.2;
    options.retinaGeodetics (1,1) logical = false
    options.rayPath {mustBeNumeric} = [];
    options.rayColor = 'red';
    options.outputRay {mustBeNumeric} = [];
    options.outputRayColor = 'red';
    options.outputRayScale (1,1) {mustBeNumeric} = 3;
    options.addLighting (1,1) logical = true
    options.viewAngle (1,2) {mustBeNumeric} = [40 40];
end



% Open a figure
if options.newFigure
    if options.visible
        figHandle = figure('Visible', 'on');
    else
        figHandle = figure('Visible', 'off');
    end
else
    figHandle = gcf;
end
hold on

% Create an empty cell array to hold the plot handles
plotHandles = gobjects(0);

% Plot the surfaceSet if provided
if ~isempty(opticalSystem)
    
    % Parse the opticalSystem
    [opticalSystem,surfaceLabels,surfaceColors] = ...
        parseOpticalSystemArgument(opticalSystem,options.surfaceSet);
    
    % Strip any nan rows from the optical system.
    opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);
    
    % Determine the number of surfaces
    nSurfaces = size(opticalSystem,1);
    
    % If the surfaceLabels or surfaceColors are unknown, fill them in
    % with unknown gray
    if isempty(surfaceColors)
        surfaceColors=cell(nSurfaces,1);
        surfaceColors(:) = {[0.5 0.5 0.5]};
    end
    if isempty(surfaceLabels)
        surfaceLabels=cell(nSurfaces,1);
        surfaceLabels(:) = {'unknown'};
    end
        
    % Loop over the surfaces
    for ii=1:nSurfaces
        % Obtain the quadric, proceed if there are no nans
        S = opticalSystem(ii,1:10);
        if ~any(isnan(S))
            % Obtain the bounding box
            boundingBox = opticalSystem(ii,12:17);
            % Obtain the surface color and set the surface alpha
            surfaceColor = surfaceColors{ii};
            % If there is a 4th value in the surface color, use this to
            % over-ride the overall surface alpha
            if size(surfaceColor,2) == 4
                surfaceAlpha = surfaceColor(4) .* options.surfaceAlpha;
                surfaceColor = surfaceColor(1:3);
            else
                surfaceAlpha = 0.5 .* options.surfaceAlpha;
            end
            % Plot the surface. If it is the retinal surface, and geodetic
            % lines have been requested, include these.
            if strcmp(surfaceLabels{ii},'retina') && options.retinaGeodetics
                plotHandles(end+1) = quadric.plotGridSurface(S,boundingBox,surfaceColor,surfaceAlpha,'g','b',options.surfaceAlpha);
            else
                plotHandles(end+1) = quadric.plotGridSurface(S,boundingBox,surfaceColor,surfaceAlpha);
            end
        end
    end
end

% Plot the rayPath if provided
if ~isempty(options.rayPath)
    plotHandles(end+1) = plot3(options.rayPath(1,:),options.rayPath(2,:),options.rayPath(3,:),'-','Color',options.rayColor);
end

% Add the outputRay if provided
if ~isempty(options.outputRay)
    outputRay = options.outputRay;
    p1=outputRay(:,1);
    p2=p1+outputRay(:,2).*options.outputRayScale;
    r = [p1 p2];
    plotHandles(end+1) = plot3(r(1,:),r(2,:),r(3,:),'-','Color',options.outputRayColor);
end

% Add a lighting source if requested and if lighting is not already present
if options.addLighting
    tmpHandle = gca;
    if ~any(isgraphics(tmpHandle.Children,'Light'))
        lightHandle = camlight;
        lighting gouraud
    end
end

% Set the viewing angle
view(options.viewAngle);

end % plotOpticalSystem

