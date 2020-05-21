function [figHandle, plotObjectHandles, renderedFrame] = renderEyePose(eyePose, sceneGeometry, varargin)
% Creates an image of the eye for a given eyePose and sceneGeometry
%
% Syntax:
%  [figHandle, plotObjectHandles, renderedFrame] = renderEyePose(eyePose, sceneGeometry)
%
% Description:
%   Create an image of the model eye given the sceneGeometry and a specific
%   eye pose.
%
% Inputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius in mm.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%   cameraTrans           - A 3x1 vector with values for [horizontal;
%                           vertical; depth] in units of mm. The values are
%                           relative to:
%                               sceneGeometry.cameraPosition.translation
%                           This key-value provides an easy way to update
%                           the camera translation.
%  'backgroundImage'      - sizeX x sizeY matrix that contains an image to
%                           displayed behind the rendered eye. Typically,
%                           this is an image of an eye for which this eye
%                           pose and sceneGeometry have been estimated.
%  'newFigure'            - Logical. Determine if we create a new figure.
%  'visible'              - Logical. If a new figure is created,
%                           determines if the figure is visible or not.
%  'showPupilTextLabels'  - Logical. Determines if the pupil perimeter
%                           points are labeled with numbers. This is used,
%                           for example, to illustrate the effects of
%                           torsion.
%  'showAzimuthPlane'     - Logical. If set to true, a magenta line is
%                           added to the rendering. The line connects the
%                           center of the pupil for a -50 and +50 rotation
%                           of azimuth and zero elevation.
%  'nStopPerimPoints'     - Scalar. The number of stop perimeter points.
%  'nIrisPerimPoints'     - Scalar. The number of iris perimeter points
%  'modelEyeLabelNames'   - Cell array of character vectors. Identifies the
%                           elements of the 'pointLabels' variable returned
%                           by projectModelEye that are to be plotted.
%                           A special case is the label 'pupilEllipse',
%                           which is not an element of pointLabels, but is
%                           recognized by this routine and prompts the
%                           plotting of the ellipse fit to the pupil
%                           perimeter.
%  'modelEyePlotColors'   - Cell array. Line spec codes for each of the
%                           elements given in modelEyeLabelNames.
%  'modelEyeAlpha'        - Scalar or vector that sets the transparency of
%                           the rendered plot elements as clear (0) or
%                           opaque (1). If passed as a vector, the length
%                           should be equal to the length of the cell
%                           arrays for modelEyeLabelNames and
%                           modelEyePlotColors. If passed as a scalar, all
%                           plot elements are assigned this alpha.
%  'modelEyeSymbolSizeScaler' - Scalar. Determines if the rendered plot
%                           elements are made larger (>1) or smaller (<1)
%                           than default.
%  'fImplicitPresent'       Logical. Older versions of MATLAB do not have
%                           the 'fImplicit' function available for plotting
%                           of implicit functions. This flag is used to
%                           indicate if the fImplicit function is
%                           available. If set to false, the ezplot function
%                           is used instead. If undefined, the routine will
%                           test if the function is available. Setting the
%                           value in the calling routine saves time over
%                           multiple calls to renderEyePose. Further, the
%                           value must be set if renderEyePose is to be
%                           called from within the parpool, as the "exist"
%                           function call that is otherwise used is not
%                           allowed in parallel execution.
%
% Outputs:
%   figHandle             - Handle to a created figure.
%   plotObjectHandles     - Array of handles to the plotted object elements
%   renderedFrame         - Structure containing the rendered frame. The
%                           field 'cdata' has the dimensions (X,Y,3), and
%                           contains the RGB image. The field 'colormap' is
%                           empty.
%
% Examples:
%{
    %% Display a 2D image of the right eye
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [-30 -5 0 3];
    renderEyePose(eyePose, sceneGeometry);
%}
%{
    %% Show the effect of eye torsion
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    renderEyePose([0 0 0 3], sceneGeometry,'showPupilTextLabels',true);
    renderEyePose([0 0 45 3], sceneGeometry,'showPupilTextLabels',true);
%}
%{
    %% Demonstrate the effect of camera translation
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [0 0 0 3];
    renderEyePose(eyePose, sceneGeometry);
    % Adjust the position in the positive x and y direction
    sceneGeometry.cameraPosition.translation = sceneGeometry.cameraPosition.translation + [5; 5; 25];
    renderEyePose(eyePose, sceneGeometry);
%}
%{
    %% Demonstrate the effect of positive camera torsion
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [0 0 0 3];
    renderEyePose(eyePose, sceneGeometry,'showPupilTextLabels',true,'nStopPerimPoints',5);
    % Adjust the camera torsion and replot
    sceneGeometry.cameraPosition.torsion = sceneGeometry.cameraPosition.torsion + 45;
    renderEyePose(eyePose, sceneGeometry,'showPupilTextLabels',true,'nStopPerimPoints',5);
%}


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('eyePose',@(x)(isnumeric(x) && all(size(x)==[1 4])));
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('cameraTrans',[0;0;0],@isvector);
p.addParameter('backgroundImage',[],@isnumeric);
p.addParameter('backgroundGray',0.5,@isscalar);
p.addParameter('newFigure',true,@islogical);
p.addParameter('visible',true,@islogical);
p.addParameter('showPupilTextLabels',false,@islogical);
p.addParameter('showAzimuthPlane',false,@islogical);
p.addParameter('nStopPerimPoints',8,@isscalar);
p.addParameter('addPseudoTorsion',true,@islogical);
p.addParameter('nIrisPerimPoints',20,@isscalar);
p.addParameter('modelEyeLabelNames', {'aziRotationCenter' 'eleRotationCenter', 'retina' 'irisPerimeter' 'stopCenter' 'pupilPerimeter' 'pupilEllipse' 'cornea' 'cornealApex' 'glint_01' 'glint_02'}, @iscell);
p.addParameter('modelEyePlotColors', {'>r' '^m' '.w' 'ob' '+r' '*g' '-g' '.y' '*y' 'xr' 'xr'}, @iscell);
p.addParameter('modelEyeAlpha',1,@isnumeric);
p.addParameter('modelEyeSymbolSizeScaler',1,@isnumeric);
p.addParameter('fImplicitPresent',[],@islogical);


% parse
p.parse(eyePose, sceneGeometry, varargin{:})

% If we have not been told if the fImplicit function is present, test for
% its existence
if isempty(p.Results.fImplicitPresent)
    fImplicitPresent = (exist('fimplicit','file')==2);
else
    fImplicitPresent = p.Results.fImplicitPresent;
end

% Expand the modelEyeAlpha parameter to the full length of the model
% components to be labeled
if length(p.Results.modelEyeAlpha)==1
    modelEyeAlpha = zeros(size(p.Results.modelEyeLabelNames));
    modelEyeAlpha(:) = p.Results.modelEyeAlpha;
else
    modelEyeAlpha = p.Results.modelEyeAlpha;
end

% Grab the image size
imageSizeX = sceneGeometry.cameraIntrinsic.sensorResolution(1);
imageSizeY = sceneGeometry.cameraIntrinsic.sensorResolution(2);

% A blank frame to initialize the figure
if ~isempty(p.Results.backgroundImage)
    backgroundImage = p.Results.backgroundImage;
    if ismatrix(backgroundImage)
        if size(backgroundImage) ~= [imageSizeY imageSizeX]
            error('renderEyePose:backgroundImageSize','The passed background image does not match the camera sensor specified in sceneGeometry');
        end
    else
        if size(backgroundImage) ~= [imageSizeY imageSizeX 3]
            error('renderEyePose:backgroundImageSize','The passed background image does not match the camera sensor specified in sceneGeometry');
        end
    end
else
    backgroundImage = zeros(imageSizeY,imageSizeX)+p.Results.backgroundGray;
end

% Open a figure
if p.Results.newFigure
    if p.Results.visible
        figHandle = figure('Visible', 'on');
    else
        figHandle = figure('Visible', 'off');
    end
    imshow(backgroundImage,[], 'Border', 'tight');
    % Prepare the figure
    axis off
    axis equal
    xlim([0 imageSizeX]);
    ylim([0 imageSizeY]);
else
    figHandle = gcf;
    % Only add a background image if one was passed
    if ~isempty(p.Results.backgroundImage)
        imshow(backgroundImage,[], 'Border', 'tight');
    end
end

% Set hold to on
hold on

% Obtain the pupilProjection of the model eye to the image plane
[pupilEllipseParams, ~, imagePoints, ~, ~, ~, pointLabels] = ...
    projectModelEye(eyePose, sceneGeometry, ...
    'cameraTrans',p.Results.cameraTrans, ...
    'fullEyeModelFlag', true, 'replaceReflectedPoints', true, ...
    'nStopPerimPoints',p.Results.nStopPerimPoints, ...
    'addPseudoTorsion',p.Results.addPseudoTorsion,...
    'nIrisPerimPoints',p.Results.nIrisPerimPoints);

% Set up an empty variable to hold plot object handles
plotObjectHandles = gobjects(0);

% Loop through the point labels present in the eye model
for pp = 1:length(p.Results.modelEyeLabelNames)
    % Check if we should plot the pupilEllipse
    if strcmp(p.Results.modelEyeLabelNames{pp},'pupilEllipse')
        % Add the pupil fit ellipse
        plotObjectHandles(end+1) = addTransparentEllipseToFigure(...
            pupilEllipseParams,imageSizeX,imageSizeY, ...
            p.Results.modelEyePlotColors{pp}(2),1,gca,fImplicitPresent);
    else
        % Plot this label
        idx = strcmp(pointLabels,p.Results.modelEyeLabelNames{pp});
        mc =  p.Results.modelEyePlotColors{pp};
        switch mc(1)
            case '.'
                s = scatter(imagePoints(idx,1), imagePoints(idx,2), (imageSizeX/(200/p.Results.modelEyeSymbolSizeScaler))^2, 'o', 'filled', 'MarkerFaceColor', mc(2), 'MarkerEdgeColor','none');
                s.MarkerFaceAlpha = modelEyeAlpha(pp);
                plotObjectHandles(end+1) = s;
            case 'o'
                s = scatter(imagePoints(idx,1), imagePoints(idx,2), (imageSizeX/(150/p.Results.modelEyeSymbolSizeScaler))^2, 'o', 'filled', 'MarkerFaceColor', mc(2), 'MarkerEdgeColor','none');
                s.MarkerFaceAlpha = modelEyeAlpha(pp);
                plotObjectHandles(end+1) = s;
            case 'O'
                s = scatter(imagePoints(idx,1), imagePoints(idx,2), (imageSizeX/(75/p.Results.modelEyeSymbolSizeScaler))^2, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', mc(2));
                s.MarkerFaceAlpha = modelEyeAlpha(pp);
                plotObjectHandles(end+1) = s;
            otherwise
                s = scatter(imagePoints(idx,1), imagePoints(idx,2), (imageSizeX/(100/p.Results.modelEyeSymbolSizeScaler))^2, mc(1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', mc(2));
                s.MarkerEdgeAlpha = modelEyeAlpha(pp);
                plotObjectHandles(end+1) = s;
        end
        % If we are plotting the pupil perimeter points, see if we would
        % like to label them
        if p.Results.showPupilTextLabels
            if strcmp(p.Results.modelEyeLabelNames{pp},'pupilPerimeterFront') || strcmp(p.Results.modelEyeLabelNames{pp},'pupilPerimeter')
                % Put text labels for the pupil perimeter points so that we can
                % follow them through rotations and translations to validate
                % the projection model
                txtHandles = text(imagePoints(idx,1), imagePoints(idx,2), num2str(find(idx)));
                plotObjectHandles = [plotObjectHandles txtHandles'];
            end
        end
        
    end
end % loop over label names

% Add a line that shows the azimuthal plane of rotation. This reflects
% camera torsion
if p.Results.showAzimuthPlane
    [~, ~, imagePoints] = projectModelEye([-50 0 0 1], sceneGeometry);
    A = nanmean(imagePoints);
    [~, ~, imagePoints] = projectModelEye([50 0 0 1], sceneGeometry);
    B = nanmean(imagePoints);
    plot([A(1) B(1)],[A(2) B(2)],'-r')
end

hold off

% Get the rendered frame
renderedFrame=getframe(figHandle);

end
