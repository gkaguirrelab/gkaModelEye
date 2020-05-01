function figHandle = plotOpticalSystem(varargin)
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
%  'surfaceAlpha'         - Scalar. The alpha transparency to use for the
%                           surfaces.
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
%
% Examples:
%{
    % Plot an entire surface set
    sceneGeometry = createSceneGeometry();
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
%}
%{
    % Plot just the opticalSystem with default color and transparency
    sceneGeometry = createSceneGeometry();
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera.opticalSystem,'addLighting',true);
%}
%{
    %% Rays from the retina through the eye and a spectacle lens
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',-2,'calcLandmarkFovea',true);
    % Plot the optical system
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    % Define an initial ray arising at the fovea
    p = sceneGeometry.eye.landmarks.fovea.coords';
    % Loop over horizontal angles relative to the visual axis
    for ii = -2:1:2
        % Assemble the ray
        R = quadric.normalizeRay( ...
            quadric.anglesToRay(p, ...
                sceneGeometry.eye.landmarks.fovea.degField(1)+ii, ...
                sceneGeometry.eye.landmarks.fovea.degField(2)));
        % Perform the ray trace
        [outputRay, rayPath] = rayTraceQuadrics(R, sceneGeometry.refraction.retinaToCamera.opticalSystem);
        % Add this ray to the optical system plot
        plotOpticalSystem('newFigure',false,'outputRay',outputRay,'rayPath',rayPath);
    end
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('newFigure',true,@islogical);
p.addParameter('visible',true,@islogical);
p.addParameter('surfaceSet',[], @(x)(isempty(x) | isstruct(x) | isnumeric(x)));
p.addParameter('surfaceAlpha', 0.1,@isnumeric);
p.addParameter('retinaGeodetics', false,@islogical);
p.addParameter('rayPath',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('rayColor','red',@(x)(ischar(x) | isnumeric(x)));
p.addParameter('outputRay',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('outputRayColor','red',@(x)(ischar(x) | isnumeric(x)));
p.addParameter('outputRayScale',3,@isnumeric);
p.addParameter('addLighting',false, @islogical);
p.addParameter('viewAngle',[40 40],@isnumeric);

% parse
p.parse(varargin{:})


% Open a figure
if p.Results.newFigure
    if p.Results.visible
        figHandle = figure('Visible', 'on');
    else
        figHandle = figure('Visible', 'off');
    end
else
    figHandle = gcf;
end
hold on

% Plot the surfaceSet if provided
if ~isempty(p.Results.surfaceSet)

    % Assemble surfaceSet components
    if isstruct(p.Results.surfaceSet)
        % If we have a whole surface set, extract the optical system and
        % plot information
        opticalSystem=p.Results.surfaceSet.opticalSystem;
        
        % Strip any nan rows from the optical system.
        opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);
        
        % Determine the number of surfaces
        nSurfaces = size(opticalSystem,1);
        
        % Obtain the surface Labels and colors
        surfaceColors=p.Results.surfaceSet.surfaceColors;
        surfaceLabels=p.Results.surfaceSet.surfaceLabels;
        
    else
        % We just have an optical system matrix.
        opticalSystem=p.Results.surfaceSet;

        % Strip any nan rows from the optical system.
        opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

        % Determine the number of surfaces
        nSurfaces = size(opticalSystem,1);

        % The surfaces are gray and the labels are unknown
        surfaceColors=cell(nSurfaces,1);
        surfaceColors(:) = {[0.5 0.5 0.5]};
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
            % Plot the surface. If it is the retinal surface, and geodetic
            % lines have been requested, include these.
            if strcmp(surfaceLabels{ii},'retina') && p.Results.retinaGeodetics
                quadric.plotSurface(S,boundingBox,surfaceColors{ii},p.Results.surfaceAlpha,'g','b',p.Results.surfaceAlpha);
            else
                quadric.plotSurface(S,boundingBox,surfaceColors{ii},p.Results.surfaceAlpha);
            end
        end
    end
end

% Plot the rayPath if provided
if ~isempty(p.Results.rayPath)
    plot3(p.Results.rayPath(1,:),p.Results.rayPath(2,:),p.Results.rayPath(3,:),'-','Color',p.Results.rayColor);
end

% Add the outputRay if provided
if ~isempty(p.Results.outputRay)
    outputRay = p.Results.outputRay;
    p1=outputRay(:,1);
    p2=p1+outputRay(:,2).*p.Results.outputRayScale;
    r = [p1 p2];
    plot3(r(1,:),r(2,:),r(3,:),'-','Color',p.Results.outputRayColor);
end

% Add a lighting source if requested
if p.Results.addLighting
    camlight
    lighting gouraud
end

% Set the viewing angle
view(p.Results.viewAngle);

end % plotOpticalSystem

