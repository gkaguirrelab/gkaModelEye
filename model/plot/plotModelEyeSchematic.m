function figHandle = plotModelEyeSchematic(eye, varargin)
% Creates a cross-section schematic illustration of the model eye
%
% Syntax:
%  figHandle = plotModelEyeSchematic(eye)
%
% Description:
%   Create a schematic diagram of the model eye specified in the passed
%   eye variable. The eye structure is described in the routine:
%       model/modelEyeParameters.m
%
% Inputs:
%   eye                   - Struct. This is the eye struct returned from
%                           modelEyeParameters(). Optionally, if a
%                           sceneGeometry structure is passed, the routine
%                           will look for the eye field.
%
% Optional key/value pairs:
%  'view'                 - String. The view to display. Valid choices
%                           include {'horizontal','vertical'};
%  'newFigure'            - Logical. Determines if we create a new figure.
%  'plotColor'            - String. Matlab line spec code for line color,
%                           e.g., {'k','r','b'};
%  'rayPath'              - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(:,1)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan. This rayPath could be obtained
%                           using the function:
%                               quadric/rayTraceQuadrics.m
%  'outputRay'            - 3x2 matrix that specifies the ray as a unit
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%  'plotIris','plotCornealApex','plotStop','plotRotationCenters',
%  'plotVisualAxis'       - Logical. Control display of these model
%                           elements.
%
% Outputs:
%   figHandle             - Handle to a created figure. Empty if a new
%                           figure was not requested.
%
% Examples:
%{
    % Basic call for a horizontal view plot
    eye = modelEyeParameters();
    plotModelEyeSchematic(eye);
%}
%{
    % Can also pass the whole sceneGeometry
    sceneGeometry = createSceneGeometry();
    plotModelEyeSchematic(sceneGeometry);
%}
%{
    % A plot with the fovea, visual axis, and line of sight
    sceneGeometry = createSceneGeometry();
    % Calculate the line-of-sight axis of the eye
    rayDestination = sceneGeometry.eye.landmarks.fovea.coords;
    rayPathLoS = calcSightRayToRetina(sceneGeometry.eye,rayDestination);
    plotModelEyeSchematic(sceneGeometry.eye,'plotVisualAxis',true,'rayPath',rayPathLoS);
%}
%{
    % Two panel plot with horizontal and vertical views for eyes with 0 and -10
    % spherical ametropia
    figure
    subplot(2,1,1)
    eye = modelEyeParameters('sphericalAmetropia',0);
    plotModelEyeSchematic(eye,'view','horizontal','newFigure',false,'plotColor','k')
    eye = modelEyeParameters('sphericalAmetropia',-10);
    plotModelEyeSchematic(eye,'view','horizontal','newFigure',false,'plotColor','r')
    subplot(2,1,2)
    eye = modelEyeParameters('sphericalAmetropia',0);
    plotModelEyeSchematic(eye,'view','vertical','newFigure',false,'plotColor','k')
    eye = modelEyeParameters('sphericalAmetropia',-10);
    plotModelEyeSchematic(eye,'view','vertical','newFigure',false,'plotColor','r')
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eye',@isstruct);

% Optional
p.addParameter('view','horizontal',@ischar);
p.addParameter('newFigure',true,@islogical);
p.addParameter('plotColor','k',@ischar);
p.addParameter('rayPath',[],@(x)(isempty(x) || ismatrix(x)));
p.addParameter('outputRay',[],@(x)(isempty(x) || ismatrix(x)));
p.addParameter('plotIris',false,@islogical);
p.addParameter('plotCornealApex',false,@islogical);
p.addParameter('plotStop',false,@islogical);
p.addParameter('plotRotationCenters',false,@islogical);
p.addParameter('plotVisualAxis',false,@islogical);
p.addParameter('plotFovea',false,@islogical);
p.addParameter('plotOpticDisc',false,@islogical);

% parse
p.parse(eye, varargin{:})

% Check if there is an eye field in the structure, in which case we are
% dealing with a sceneGeometry, and we will want to extract the eye field.
if isfield(eye,'eye')
    eye = eye.eye;
end

% Open a figure
if p.Results.newFigure
    figHandle = figure;
end
hold on

% Determine which view we are using
switch p.Results.view
    case {'axial','Axial','Ax','ax','Horizontal','horizontal','horiz'}
        PdimA = 1;
        PdimB = 2;
        titleString = 'Horizontal';
        switch eye.meta.eyeLaterality
            case 'Left'
                yLabelString = 'temporal <----> nasal';
            case 'Right'
                yLabelString = 'nasal <----> temporal';
        end
        xLabelString = 'posterior <----> anterior';
    case {'sagittal','Sagittal','Sag','sag','Vertical','vertical','vert'}
        PdimA = 1;
        PdimB = 3;
        titleString = 'Vertical';
        yLabelString = 'inferior <----> superior';
        xLabelString = 'posterior <----> anterior';
    otherwise
        error('Not a recognized view for the schematic eye');
end

%% Plot the anterior chamber, vitreous chamber, and the lens surfaces
quadric.plotConicSection(eye.retina.S(1:10), titleString, p.Results.plotColor, eye.retina.boundingBox)
quadric.plotConicSection(eye.cornea.S(1,:), titleString, p.Results.plotColor, eye.cornea.boundingBox(1,:))
quadric.plotConicSection(eye.cornea.S(end,:), titleString, p.Results.plotColor, eye.cornea.boundingBox(end,:))
quadric.plotConicSection(eye.lens.S(1,:), titleString, p.Results.plotColor, eye.lens.boundingBox(1,:))
quadric.plotConicSection(eye.lens.S(end,:), titleString, p.Results.plotColor, eye.lens.boundingBox(end,:))


%% Add various additional elements, under the control of flags
if p.Results.plotStop
    plot([eye.stop.center(PdimA) eye.stop.center(PdimA)],[-2 2],['-' p.Results.plotColor]);
end

if p.Results.plotRotationCenters
    plot(eye.rotationCenters.azi(PdimA),eye.rotationCenters.azi(PdimB),['*' p.Results.plotColor])
    plot(eye.rotationCenters.ele(PdimA),eye.rotationCenters.ele(PdimB),['o' p.Results.plotColor])
end

if p.Results.plotIris
    plot(eye.iris.center(PdimA),eye.iris.center(PdimB)+eye.iris.radius,['x' p.Results.plotColor])
    plot(eye.iris.center(PdimA),eye.iris.center(PdimB)-eye.iris.radius,['x' p.Results.plotColor])
end

if p.Results.plotCornealApex
    sg.eye = eye;
    [~, ~, ~, ~, ~, eyeWorldPoints, pointLabels] = projectModelEye([0 0 0 1], sg, 'fullEyeModelFlag',true);
    idx = find(strcmp(pointLabels,'cornealApex'));
    plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),['*' p.Results.plotColor]);
end

if isfield(eye,'landmarks')
    
    % fovea
    if isfield(eye.landmarks,'fovea')
        if p.Results.plotFovea
        plot(eye.landmarks.fovea.coords(PdimA),eye.landmarks.fovea.coords(PdimB),['*' p.Results.plotColor])
        end        
        if p.Results.plotVisualAxis
            % Obtain the nodal ray from the fovea
            rayPath = calcNodalRayToRetina(eye,eye.landmarks.fovea.coords);
            plot(rayPath(PdimA,:),rayPath(PdimB,:),['-' p.Results.plotColor]);
        end
    end
    
    % optic disc
    if isfield(eye.landmarks,'opticDisc')
        if p.Results.plotOpticDisc
        plot(eye.landmarks.opticDisc.coords(PdimA),eye.landmarks.opticDisc.coords(PdimB),['x' p.Results.plotColor])
        end
    end
end


%% Plot a passed rayPath
if ~isempty(p.Results.rayPath)
    plot(p.Results.rayPath(PdimA,:),p.Results.rayPath(PdimB,:),['-' p.Results.plotColor]);
end


%% Plot a passed outputRay
if ~isempty(p.Results.outputRay)
    p1=p.Results.outputRay(:,1);
    p2=p1+p.Results.outputRay(:,2).*3;
    r = [p1 p2];
    plot(r(PdimA,:),r(PdimB,:),['-' p.Results.plotColor]);
end


%% Reference axis
xRange = xlim;
plot(xRange,[0 0],['-' p.Results.plotColor]);

axis equal
title(titleString);
ylabel(yLabelString);
xlabel(xLabelString);

%% Set the xlim
xlim([eye.landmarks.vertex.coords(1)-5 5])
ylim([-abs(eye.landmarks.vertex.coords(1)/2) abs(eye.landmarks.vertex.coords(1)/2)])
    
end

