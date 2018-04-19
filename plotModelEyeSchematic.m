function figHandle = plotModelEyeSchematic(eye, varargin)
% Creates a cross-section schematic illustration of the model eye
%
% Syntax:
%  plotModelEyeSchematic(eye)
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
    eye = modelEyeParameters;
    plotModelEyeSchematic(eye)
%}
%{
    % Two panel plot with axial and sagittal views for eyes with 0 and -10
    % spherical ametropia
    figure
    subplot(2,1,1)
    eye = modelEyeParameters('sphericalAmetropia',0);
    plotModelEyeSchematic(eye,'view','axial','newFigure',false,'plotColor','k')
    eye = modelEyeParameters('sphericalAmetropia',-10);
    plotModelEyeSchematic(eye,'view','axial','newFigure',false,'plotColor','r')
    subplot(2,1,2)
    eye = modelEyeParameters('sphericalAmetropia',0);
    plotModelEyeSchematic(eye,'view','sagittal','newFigure',false,'plotColor','k')
    eye = modelEyeParameters('sphericalAmetropia',-10);
    plotModelEyeSchematic(eye,'view','sagittal','newFigure',false,'plotColor','r')
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
        SdimA = 3;
        SdimB = 1;
        rotationField = 'azi';
        titleString = 'Axial';
        yLabelString = 'temporal <----> nasal';
        xLabelString = 'posterior <----> anterior';
        postChamberRange = [-30, -7, -15, 15];
        corneaRange = [eye.pupil.center(1), 5, -15, 15];
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
        corneaRange = [eye.pupil.center(1), 5, -15, 15];
        lensRange = [-10, 0, -5, 5];       
    otherwise
        error('Not a recognized view for the schematic eye');
end

%% Plot the anterior and posterior chambers and the lens
ep = ellipse_ex2im([eye.posteriorChamber.center([PdimA PdimB]) eye.posteriorChamber.radii([PdimA PdimB]) 0]);
plotEllipse(ep,p.Results.plotColor,postChamberRange)
ep = ellipse_ex2im([eye.cornea.back.center([PdimA PdimB]) eye.cornea.back.radii([PdimA PdimB]) deg2rad(eye.cornea.axis(PdimA))]);
plotEllipse(ep,p.Results.plotColor,corneaRange)
ep = ellipse_ex2im([eye.cornea.front.center([PdimA PdimB]) eye.cornea.front.radii([PdimA PdimB]) deg2rad(eye.cornea.axis(PdimA))]);
plotEllipse(ep,p.Results.plotColor,corneaRange)
ep = ellipse_ex2im([eye.lens.front.center([PdimA PdimB]) eye.lens.front.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,p.Results.plotColor,lensRange)
ep = ellipse_ex2im([eye.lens.back.center([PdimA PdimB]) eye.lens.back.radii([PdimA PdimB]) 0]);
plotHyperbola(ep,p.Results.plotColor,lensRange)

%% Add a 2mm radius pupil, center of rotation, iris boundary, and fovea
plot([eye.pupil.center(PdimA) eye.pupil.center(PdimA)],[-2 2],['-' p.Results.plotColor]);
plot(eye.rotationCenters.azi(PdimA),eye.rotationCenters.azi(PdimB),['>' p.Results.plotColor])
plot(eye.rotationCenters.ele(PdimA),eye.rotationCenters.ele(PdimB),['^' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)+eye.iris.radius,['x' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)-eye.iris.radius,['x' p.Results.plotColor])
plot(eye.posteriorChamber.fovea(PdimA),eye.posteriorChamber.fovea(PdimB),['*' p.Results.plotColor])

%% Plot the cornealApex
sg.eye = eye;
[~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd([0 0 0 1], sg, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'cornealApex'));
plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),['*' p.Results.plotColor]);

%% Plot the visual axis
m = (eye.posteriorChamber.fovea(PdimB) - eye.lens.nodalPoint.rear(PdimB)) / (eye.posteriorChamber.fovea(PdimA) - eye.lens.nodalPoint.rear(PdimA));
b = eye.lens.nodalPoint.rear(PdimB) -  (eye.lens.nodalPoint.rear(PdimA) * m);
xRange = xlim;
plot(xRange,xRange.*m+b,[':' p.Results.plotColor]);


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
