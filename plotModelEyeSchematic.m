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
%  'crossSectionView'     - String. The view to display. Valid choices
%                           include {'axial','sagittal','coronal'};
%  'newFigure'            - Logical. Determines if we create a new figure.
%
% Outputs:
%   figHandle             - Handle to a created figure. Empty if a new
%                           figure was not requested.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry;
    plotModelEyeSchematic(sceneGeometry)
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eye',@isstruct);

% Optional
p.addParameter('crossSectionView','axial',@ischar);
p.addParameter('newFigure',true,@islogical);

% parse
p.parse(sceneGeometry, varargin{:})


% Open a figure
if p.Results.newFigure
    figHandle = figure;
end
hold on

% Determine which view we are using
switch p.Results.crossSectionView
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
    case 'coronal'
        PdimA = 2;
        PdimB = 3;
        rotationField = 'tor';
        titleString = 'Coronal';
        xLabelString = 'temporal <----> nasal';
        yLabelString = 'inferior <----> superior';
    otherwise
        error('Not a recognized crossSectionView for the schematic eye');
end

%% Plot the anterior and posterior chambers and the lens
ep = ellipse_ex2im([sceneGeometry.eye.posteriorChamber.center([PdimA PdimB]) sceneGeometry.eye.posteriorChamber.radii([PdimA PdimB]) 0]);
plotEllipse(ep,'k',postChamberRange)
ep = ellipse_ex2im([sceneGeometry.eye.cornea.back.center([PdimA PdimB]) sceneGeometry.eye.cornea.back.radii([PdimA PdimB]) deg2rad(sceneGeometry.eye.cornea.axis(PdimA))]);
plotEllipse(ep,'k',corneaRange)
ep = ellipse_ex2im([sceneGeometry.eye.cornea.front.center([PdimA PdimB]) sceneGeometry.eye.cornea.front.radii([PdimA PdimB]) deg2rad(sceneGeometry.eye.cornea.axis(PdimA))]);
plotEllipse(ep,'k',corneaRange)
ep = ellipse_ex2im([sceneGeometry.eye.lens.front.center([PdimA PdimB]) sceneGeometry.eye.lens.front.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,[0.5 0.5 0.5],lensRange)
ep = ellipse_ex2im([sceneGeometry.eye.lens.back.center([PdimA PdimB]) sceneGeometry.eye.lens.back.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,[0.5 0.5 0.5],lensRange)

%% Add the pupil center, 2mm radius pupil, center of rotation, iris boundary
plot([sceneGeometry.eye.pupil.center(PdimA) sceneGeometry.eye.pupil.center(PdimA)],[-2 2],'-k');
plot(sceneGeometry.eye.pupil.center(PdimA),sceneGeometry.eye.pupil.center(PdimB),'*r')
plot(sceneGeometry.eye.rotationCenters.azi(PdimA),sceneGeometry.eye.rotationCenters.azi(PdimB),'>g')
plot(sceneGeometry.eye.rotationCenters.ele(PdimA),sceneGeometry.eye.rotationCenters.ele(PdimB),'^g')
plot(sceneGeometry.eye.iris.center(PdimA),sceneGeometry.eye.iris.center(PdimB)+sceneGeometry.eye.iris.radius,'xb')
plot(sceneGeometry.eye.iris.center(PdimA),sceneGeometry.eye.iris.center(PdimB)-sceneGeometry.eye.iris.radius,'xb')

%% Plot the cornealApex
[~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd([0 0 0 1], sceneGeometry, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'cornealApex'));
plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),'xy');

%% Plot the fixation Axis
[~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd([sceneGeometry.eye.gamma(1) sceneGeometry.eye.gamma(2) sceneGeometry.eye.gamma(3) 1], sceneGeometry, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'opticalAxisOrigin'));
rotationCenter = mean([sceneGeometry.eye.rotationCenters.azi; sceneGeometry.eye.rotationCenters.ele]);
plot([sceneWorldPoints(idx,SdimA) rotationCenter(PdimA)],[sceneWorldPoints(idx,SdimB) rotationCenter(PdimB)],'-r');

%% Reference axis
refline(0,0)

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
