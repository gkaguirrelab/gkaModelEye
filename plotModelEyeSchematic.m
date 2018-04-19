function figHandle = plotModelEyeSchematic(sceneGeometry, varargin)
% Creates a cross-section schematic illustration of the model eye
%
% Syntax:
%  plotModelEyeSchematic(eye)
%
% Description:
%   Create a schematic diagram of the model eye specified in the passed
%   sceneGeometry variable.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry struct.
%
% Optional key/value pairs:
%  'view'     - String. The view to display. Valid choices
%                           include {'axial','sagittal'};
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
    % Basic call for an axial view plot in a new window with black lines
    sceneGeometry = createSceneGeometry;
    plotModelEyeSchematic(sceneGeometry)
%}
%{
    % Two panel plot with axial and sagittal views for eyes with 0 and -10
    % spherical ametropia
    figure
    subplot(2,1,1)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',0);
    plotModelEyeSchematic(sceneGeometry,'view','axial','newFigure',false,'plotColor','k')
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-10);
    plotModelEyeSchematic(sceneGeometry,'view','axial','newFigure',false,'plotColor','r')
    subplot(2,1,2)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',0);
    plotModelEyeSchematic(sceneGeometry,'view','sagittal','newFigure',false,'plotColor','k')
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-10);
    plotModelEyeSchematic(sceneGeometry,'view','sagittal','newFigure',false,'plotColor','r')
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eye',@isstruct);

% Optional
p.addParameter('view','axial',@ischar);
p.addParameter('newFigure',true,@islogical);
p.addParameter('plotColor','k',@ischar);

% parse
p.parse(sceneGeometry, varargin{:})


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
        SdimA = 3;
        SdimB = 1;
        rotationField = 'azi';
        titleString = 'Axial';
        yLabelString = 'temporal <----> nasal';
        xLabelString = 'posterior <----> anterior';
        postChamberRange = [-30, -7, -15, 15];
        corneaRange = [sceneGeometry.eye.pupil.center(1), 5, -15, 15];
        lensRange = [-10, 0, -5, 5];
    case {'sagittal','Sagittal','Sag','sag','Vertical','vertical','vert'}
        PdimA = 1;
        PdimB = 3;
        SdimA = 3;
        SdimB = 2;
        set(gca,'Ydir','reverse')
        rotationField = 'ele';
        titleString = 'Sagittal';
        yLabelString = 'inferior <----> superior';
        xLabelString = 'posterior <----> anterior';
        postChamberRange = [-30, -7, -15, 15];
        corneaRange = [sceneGeometry.eye.pupil.center(1), 5, -15, 15];
        lensRange = [-10, 0, -5, 5];       
    otherwise
        error('Not a recognized view for the schematic eye');
end

%% Plot the anterior and posterior chambers and the lens
ep = ellipse_ex2im([sceneGeometry.eye.posteriorChamber.center([PdimA PdimB]) sceneGeometry.eye.posteriorChamber.radii([PdimA PdimB]) 0]);
plotEllipse(ep,p.Results.plotColor,postChamberRange)
ep = ellipse_ex2im([sceneGeometry.eye.cornea.back.center([PdimA PdimB]) sceneGeometry.eye.cornea.back.radii([PdimA PdimB]) deg2rad(sceneGeometry.eye.cornea.axis(PdimA))]);
plotEllipse(ep,p.Results.plotColor,corneaRange)
ep = ellipse_ex2im([sceneGeometry.eye.cornea.front.center([PdimA PdimB]) sceneGeometry.eye.cornea.front.radii([PdimA PdimB]) deg2rad(sceneGeometry.eye.cornea.axis(PdimA))]);
plotEllipse(ep,p.Results.plotColor,corneaRange)
ep = ellipse_ex2im([sceneGeometry.eye.lens.front.center([PdimA PdimB]) sceneGeometry.eye.lens.front.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,p.Results.plotColor,lensRange)
ep = ellipse_ex2im([sceneGeometry.eye.lens.back.center([PdimA PdimB]) sceneGeometry.eye.lens.back.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,p.Results.plotColor,lensRange)

%% Add a 2mm radius pupil, center of rotation, iris boundary, and fovea
plot([sceneGeometry.eye.pupil.center(PdimA) sceneGeometry.eye.pupil.center(PdimA)],[-2 2],['-' p.Results.plotColor]);
plot(sceneGeometry.eye.rotationCenters.azi(PdimA),sceneGeometry.eye.rotationCenters.azi(PdimB),['>' p.Results.plotColor])
plot(sceneGeometry.eye.rotationCenters.ele(PdimA),sceneGeometry.eye.rotationCenters.ele(PdimB),['^' p.Results.plotColor])
plot(sceneGeometry.eye.iris.center(PdimA),sceneGeometry.eye.iris.center(PdimB)+sceneGeometry.eye.iris.radius,['x' p.Results.plotColor])
plot(sceneGeometry.eye.iris.center(PdimA),sceneGeometry.eye.iris.center(PdimB)-sceneGeometry.eye.iris.radius,['x' p.Results.plotColor])
plot(sceneGeometry.eye.posteriorChamber.fovea(PdimA),sceneGeometry.eye.posteriorChamber.fovea(PdimB),['*' p.Results.plotColor])

%% Plot the cornealApex
[~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd([0 0 0 1], sceneGeometry, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'cornealApex'));
plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),['*' p.Results.plotColor]);

%% Plot the visual axis
m = (sceneGeometry.eye.posteriorChamber.fovea(PdimB) - sceneGeometry.eye.lens.nodalPoint.rear(PdimB)) / (sceneGeometry.eye.posteriorChamber.fovea(PdimA) - sceneGeometry.eye.lens.nodalPoint.rear(PdimA));
b = sceneGeometry.eye.lens.nodalPoint.rear(PdimB) -  (sceneGeometry.eye.lens.nodalPoint.rear(PdimA) * m);
xRange = xlim;
plot(xRange,xRange.*m+b,[':' p.Results.plotColor]);

%% Plot the visual axis when the eye is rotated to -alpha
% This can be used to confirm that we are calculating alpha properly, but
% usually does not need to be displayed
%{
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd([-sceneGeometry.eye.alpha(1) -sceneGeometry.eye.alpha(2) -sceneGeometry.eye.alpha(3) 1], sceneGeometry, 'fullEyeModelFlag',true);
    idx1 = find(strcmp(pointLabels,'fovea'));
    idx2 = find(strcmp(pointLabels,'nodalPointRear'));
    m = (sceneWorldPoints(idx2,SdimB) - sceneWorldPoints(idx1,SdimB)) / (sceneWorldPoints(idx2,SdimA) - sceneWorldPoints(idx1,SdimA));
    b = sceneWorldPoints(idx1,SdimB) -  (sceneWorldPoints(idx1,SdimA) * m);
    xRange = xlim;
    plot(xRange,xRange.*m+b,[':' p.Results.plotColor]);
%}

%% Reference axis
plot(xRange,[0 0],['-' p.Results.plotColor]);

axis equal
title(titleString);
ylabel(yLabelString);
xlabel(xLabelString);

end

function plotEllipse(ep,colorCode,rangeVec)
fh=@(x,y) ep(1).*x.^2 +ep(2).*x.*y +ep(3).*y.^2 +ep(4).*x +ep(5).*y +ep(6);
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end

function plotHyperbola(ep,colorCode,rangeVec)
fh=@(x,y) ep(1).*x.^2 +ep(2).*x.*y -ep(3).*y.^2 +ep(4).*x -ep(5).*y +ep(6);
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end
