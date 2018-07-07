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
%  'view'                 - String. The view to display. Valid choices
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
    % Basic call for an axial and sagittal view plot in new windows
    eye = modelEyeParameters;
    plotModelEyeSchematic(eye);
    plotModelEyeSchematic(eye,'view','sag');
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
    case {'sagittal','Sagittal','Sag','sag','Vertical','vertical','vert'}
        PdimA = 1;
        PdimB = 3;
        SdimA = 3;
        SdimB = 2;
        rotationField = 'ele';
        titleString = 'Sagittal';
        yLabelString = 'inferior <----> superior';
        xLabelString = 'posterior <----> anterior';
    otherwise
        error('Not a recognized view for the schematic eye');
end

%% Plot the anterior and posterior chambers and the lens
plotConicSection(eye.posteriorChamber.S, titleString, p.Results.plotColor,eye.posteriorChamber.boundingBox)
plotConicSection(eye.cornea.back.S, titleString, p.Results.plotColor,eye.cornea.back.boundingBox)
plotConicSection(eye.cornea.front.S, titleString, p.Results.plotColor,eye.cornea.front.boundingBox)
plotConicSection(eye.lens.front.S, titleString, p.Results.plotColor,eye.lens.front.boundingBox)
plotConicSection(eye.lens.back.S, titleString, p.Results.plotColor,eye.lens.back.boundingBox)


%% Add a 2mm radius pupil, center of rotation, iris boundary, fovea, and optic disc
plot([eye.pupil.center(PdimA) eye.pupil.center(PdimA)],[-2 2],['-' p.Results.plotColor]);
plot(eye.rotationCenters.azi(PdimA),eye.rotationCenters.azi(PdimB),['>' p.Results.plotColor])
plot(eye.rotationCenters.ele(PdimA),eye.rotationCenters.ele(PdimB),['^' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)+eye.iris.radius,['x' p.Results.plotColor])
plot(eye.iris.center(PdimA),eye.iris.center(PdimB)-eye.iris.radius,['x' p.Results.plotColor])
plot(eye.posteriorChamber.fovea(PdimA),eye.posteriorChamber.fovea(PdimB),['*' p.Results.plotColor])
plot(eye.posteriorChamber.opticDisc(PdimA),eye.posteriorChamber.opticDisc(PdimB),['x' p.Results.plotColor])

%% Plot the cornealApex
sg.eye = eye;
[~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd([0 0 0 1], sg, 'fullEyeModelFlag',true);
idx = find(strcmp(pointLabels,'cornealApex'));
plot(eyeWorldPoints(idx,PdimA),eyeWorldPoints(idx,PdimB),['*' p.Results.plotColor]);

%% Plot the visual axis
m = (eye.posteriorChamber.fovea(PdimB) - eye.lens.nodalPoint(PdimB)) / (eye.posteriorChamber.fovea(PdimA) - eye.lens.nodalPoint(PdimA));
b = eye.lens.nodalPoint(PdimB) -  (eye.lens.nodalPoint(PdimA) * m);
xRange = xlim;
plot(xRange,xRange.*m+b,[':' p.Results.plotColor]);

%% Plot the blind spot axis
m = (eye.posteriorChamber.opticDisc(PdimB) - eye.lens.nodalPoint(PdimB)) / (eye.posteriorChamber.opticDisc(PdimA) - eye.lens.nodalPoint(PdimA));
b = eye.lens.nodalPoint(PdimB) -  (eye.lens.nodalPoint(PdimA) * m);
xRange = xlim;
plot(xRange,xRange.*m+b,[':' p.Results.plotColor]);


%% Reference axis
plot(xRange,[0 0],['-' p.Results.plotColor]);

axis equal
title(titleString);
ylabel(yLabelString);
xlabel(xLabelString);

end

function plotConicSection(S,plane,colorCode,boundingBox)
    F = quadric.vecToFunc(S);
    switch plane
        case 'Axial'
            fh = @(x,y) F(x,y,0);
            rangeVec = boundingBox([1 2 3 4]);
        case 'Sagittal'
            fh = @(x,y) F(x,0,y);
            rangeVec = boundingBox([1 2 5 6]);
        case 'Coronal'
            fh = @(x,y) F(0,x,y);
            rangeVec = boundingBox([3 4 5 6]);
    end
    fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end

function plotEllipse(ep,colorCode,rangeVec)
fh=@(x,y) ep(1).*x.^2 +ep(2).*x.*y +ep(3).*y.^2 +ep(4).*x +ep(5).*y +ep(6);
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end

function plotHyperbola(ep,colorCode,rangeVec)
fh=@(x,y) ep(1).*x.^2 +ep(2).*x.*y -ep(3).*y.^2 +ep(4).*x -ep(5).*y +ep(6);
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end
