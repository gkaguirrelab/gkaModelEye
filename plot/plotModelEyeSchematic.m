function figHandle = plotModelEyeSchematic(eye, varargin)
% Creates a cross-section schematic illustration of the model eye
%
% Syntax:
%  figHandle = plotModelEyeSchematic(eye)
%
% Description:
%   Create a schematic diagram of the model eye specified in the passed
%   eye variable.
%
% Inputs:
%   eye                   - An eye struct returned from
%                           modelEyeParameters()
%
% Optional key/value pairs:
%  'view'                 - String. The view to display. Valid choices
%                           include {'horizontal','vertical'};
%  'newFigure'            - Logical. Determines if we create a new figure.
%  'plotColor'            - String. Matlab line spec code for line color,
%                           e.g., {'k','r','b'};
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
    % A plot with the fovea, visual axis, and optical center
    eye = modelEyeParameters('calcLandmarkFovea',true,'calcLandmarkOpticalCenter',true);
    plotModelEyeSchematic(eye);
%}
%{
    % A plot with the fovea, visual axis, and line of sight
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    [outputRayLoS,rayPathLoS] = calcLineOfSightRay(sceneGeometry);
    plotModelEyeSchematic(sceneGeometry.eye,'rayPath',rayPathLoS,'outputRay',outputRayLoS);
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

% parse
p.parse(eye, varargin{:})

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
        yLabelString = 'temporal <----> nasal';
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
plotConicSection(eye.retina.S(1:10), titleString, p.Results.plotColor, eye.retina.boundingBox)
plotConicSection(eye.cornea.S(1,:), titleString, p.Results.plotColor, eye.cornea.boundingBox(1,:))
plotConicSection(eye.cornea.S(end,:), titleString, p.Results.plotColor, eye.cornea.boundingBox(end,:))
plotConicSection(eye.lens.S(1,:), titleString, p.Results.plotColor, eye.lens.boundingBox(1,:))
plotConicSection(eye.lens.S(end,:), titleString, p.Results.plotColor, eye.lens.boundingBox(end,:))


%% Add a 2mm radius pupil, center of rotation, iris boundary, fovea, and optic disc
plot([eye.stop.center(PdimA) eye.stop.center(PdimA)],[-2 2],['-' p.Results.plotColor]);
plot(eye.rotationCenters.azi(PdimA),eye.rotationCenters.azi(PdimB),['>' p.Results.plotColor])
plot(eye.rotationCenters.ele(PdimA),eye.rotationCenters.ele(PdimB),['^' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)+eye.iris.radius,['x' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)-eye.iris.radius,['x' p.Results.plotColor])

%% Plot the cornealApex
sg.eye = eye;
[~, ~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd([0 0 0 1], sg, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'cornealApex'));
plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),['*' p.Results.plotColor]);

%% Plot the fovea, the visual axis, and the line of sight
% Obtain the rayPath through the optical system from the fovea to cornea
if isfield(eye,'landmarks')
    if isfield(eye.landmarks,'fovea')
        plot(eye.landmarks.fovea.coords(PdimA),eye.landmarks.fovea.coords(PdimB),['*' p.Results.plotColor])
        
        % Obtain the nodal ray from the fovea
        [outputRay,rayPath] = calcNodalRay(eye,[],eye.landmarks.fovea.coords);
        plot(rayPath(PdimA,:),rayPath(PdimB,:),[':' p.Results.plotColor]);
        p1=outputRay(:,1);
        p2=p1+outputRay(:,2).*3;
        r = [p1 p2];
        plot(r(PdimA,:),r(PdimB,:),[':' p.Results.plotColor]);
    end
end

%% Plot the optic disc and blind spot axis
% Obtain the rayPath through the optical system from the opticDisc to cornea
if isfield(eye,'landmarks')
    if isfield(eye.landmarks,'opticDisc')
        plot(eye.landmarks.opticDisc.coords(PdimA),eye.landmarks.opticDisc.coords(PdimB),['x' p.Results.plotColor])
    end
end

%% Plot the optical center
if isfield(eye,'landmarks')
    if isfield(eye.landmarks,'opticalCenter')
        plot(eye.landmarks.opticalCenter(PdimA),eye.landmarks.opticalCenter(PdimB),['o' p.Results.plotColor]);
    end
end

%% Plot a passed rayPath
if ~isempty(p.Results.rayPath)
        plot(p.Results.rayPath(PdimA,:),p.Results.rayPath(PdimB,:),['--' p.Results.plotColor]);
end

%% Plot a passed outputRay
if ~isempty(p.Results.outputRay)
        p1=p.Results.outputRay(:,1);
        p2=p1+p.Results.outputRay(:,2).*3;
        r = [p1 p2];
        plot(r(PdimA,:),r(PdimB,:),['--' p.Results.plotColor]);
end


%% Reference axis
xRange = xlim;
plot(xRange,[0 0],['-' p.Results.plotColor]);

axis equal
title(titleString);
ylabel(yLabelString);
xlabel(xLabelString);

end

function plotConicSection(S,plane,colorCode,boundingBox)
F = quadric.vecToFunc(S);
switch plane
    case 'Horizontal'
        fh = @(x,y) F(x,y,0);
        rangeVec = boundingBox([1 2 3 4]);
    case 'Vertical'
        fh = @(x,y) F(x,0,y);
        rangeVec = boundingBox([1 2 5 6]);
    case 'Coronal'
        fh = @(x,y) F(0,x,y);
        rangeVec = boundingBox([3 4 5 6]);
end
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end
