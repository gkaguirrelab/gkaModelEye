function [figHandle, renderedFrame] = renderEyePose(eyePose, sceneGeometry, varargin)
% Creates an image of the eye for a given eyePose and sceneGeometry
%
% Syntax:
%  figHandle = renderEyePose(eyePose, sceneGeometry)
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
%  'nPupilPerimPoints'    - Scalar. The number of pupil perimeter points.
%  'nIrisPerimPoints'     - Scalar. The number of iris perimeter points
%  'modelEyeLabelNames'   - Cell array of character vectors. Identifies the
%                           elements of the 'pointLabels' variable returned
%                           by pupilProjection_fwd that are to be plotted.
%                           A special case is the label 'pupilEllipse',
%                           which is not an element of pointLabels, but is
%                           recognized by this routine and prompts the
%                           plotting of the ellipse fit to the pupil
%                           perimeter.
%  'modelEyePlotColors'   - Cell array. Line spec codes for each of the
%                           elements given in modelEyeLabelNames.
%
% Outputs:
%   figHandle             - Handle to a created figure.
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
    %% Plot the pupil ellipse for various eye poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Prepare a figure
    modelEyeLabelNames = {'pupilPerimeterBack_hidden' 'pupilPerimeterBack' 'pupilEllipse' 'pupilPerimeterFront_hidden' 'pupilPerimeterFront'};
	modelEyePlotColors = {'xr' '*r' '-y' 'xg' '*g'};    
    renderEyePose([0 0 0 2], sceneGeometry,'modelEyeLabelNames',modelEyeLabelNames,'modelEyePlotColors',modelEyePlotColors);
    for azi = -35:35:35
        for ele = -35:35:35
            eyePose = [azi ele 0 2];
            renderEyePose(eyePose, sceneGeometry,'newFigure',false,'modelEyeLabelNames',modelEyeLabelNames,'modelEyePlotColors',modelEyePlotColors);
        end
    end
%}
%{
    %% Show the effect of eye torsion
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    renderEyePose([0 0 0 3], sceneGeometry,'showPupilTextLabels',true,'nPupilPerimPoints',5);
    renderEyePose([0 0 45 3], sceneGeometry,'showPupilTextLabels',true,'nPupilPerimPoints',5);
%}
%{
    %% Demonstrate the effect of camera translation
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [0 0 0 3];
    renderEyePose(eyePose, sceneGeometry);
    % Adjust the position in the positive x and y direction
    sceneGeometry.cameraPosition.translation = sceneGeometry.cameraPosition.translation + [5; 5; 0];
    renderEyePose(eyePose, sceneGeometry);
%}
%{
    %% Demonstrate the effect of positive camera torsion
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [0 0 0 3];
    renderEyePose(eyePose, sceneGeometry,'showPupilTextLabels',true,'nPupilPerimPoints',5);
    % Adjust the camera torsion and replot
    sceneGeometry.cameraPosition.torsion = sceneGeometry.cameraPosition.torsion + 45;
    renderEyePose(eyePose, sceneGeometry,'showPupilTextLabels',true,'nPupilPerimPoints',5);
%}


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('eyePose',@(x)(isnumeric(x) && all(size(x)==[1 4])));
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('backgroundImage',[],@isnumeric);
p.addParameter('newFigure',true,@islogical);
p.addParameter('visible',true,@islogical);
p.addParameter('showPupilTextLabels',false,@islogical);
p.addParameter('nPupilPerimPoints',8,@isnumeric);
p.addParameter('nIrisPerimPoints',20,@isnumeric);
p.addParameter('modelEyeLabelNames', {'aziRotationCenter', 'eleRotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeterBack' 'pupilEllipse' 'pupilPerimeterFront' 'anteriorChamber' 'cornealApex'}, @iscell);
p.addParameter('modelEyePlotColors', {'>r' '^m' '.w' '.b' '*g' '-g' '*g' '.y' '*y'}, @iscell);
p.addParameter('modelEyeAlpha',1,@isnumeric);

% parse
p.parse(eyePose, sceneGeometry, varargin{:})

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

% A blank frame to initialize each frame
if ~isempty(p.Results.backgroundImage)
    backgroundImage = p.Results.backgroundImage;
    if size(backgroundImage) ~= [imageSizeX imageSizeY]
        error('renderEyePose:backgroundImageSize','The passed background image does not match the camera sensor specified in sceneGeometry');
    end
else
    backgroundImage = zeros(imageSizeY,imageSizeX)+0.5;
end

% Open a figure
if p.Results.newFigure
    if p.Results.visible
        figHandle = figure('Visible', 'on');
    else
        figHandle = figure('Visible', 'off');
    end
    imshow(backgroundImage, 'Border', 'tight');
    % Prepare the figure
    hold on
    axis off
    axis equal
    xlim([0 imageSizeX]);
    ylim([0 imageSizeY]);
else
    figHandle = gcf;
end

% Silence a ray tracing warning that can occur when the eye is rotated at a
% large angle and points in the iris (and sometimes pupil) encounter
% internal reflection
warnState = warning();
warning('Off','rayTraceEllipsoids:criticalAngle');

% Obtain the pupilProjection of the model eye to the image plane
[pupilEllipseParams, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd(eyePose, sceneGeometry, 'fullEyeModelFlag', true, 'nPupilPerimPoints',p.Results.nPupilPerimPoints, 'nIrisPerimPoints',p.Results.nIrisPerimPoints);

% Restore the warning state
warning(warnState);

% Loop through the point labels present in the eye model
for pp = 1:length(p.Results.modelEyeLabelNames)
    % Check if we should plot the pupilEllipse
    if strcmp(p.Results.modelEyeLabelNames{pp},'pupilEllipse')
        % Add the pupil fit ellipse
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilEllipseParams));
        fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
        % Superimpose the ellipse using fimplicit or ezplot (ezplot is the
        % fallback option for older Matlab versions)
        if exist('fimplicit','file')==2
            fimplicit(fh,[1, imageSizeX, 1, imageSizeY],'Color', p.Results.modelEyePlotColors{pp}(2),'LineWidth',1);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis off;
        else
            plotHandle=ezplot(fh,[1, imageSizeX, 1, imageSizeY]);
            set(plotHandle, 'Color', p.Results.modelEyePlotColors{pp}(2))
            set(plotHandle,'LineWidth',1);
        end
    else
        % Plot this label
        idx = strcmp(pointLabels,p.Results.modelEyeLabelNames{pp});
        mc =  p.Results.modelEyePlotColors{pp};
        switch mc(1)
            case '.'
                sc = scatter(imagePoints(idx,1), imagePoints(idx,2), 10, 'o', 'filled', 'MarkerFaceColor', mc(2), 'MarkerEdgeColor','none');
                sc.MarkerFaceAlpha = modelEyeAlpha(pp);
            case 'o'
                sc = scatter(imagePoints(idx,1), imagePoints(idx,2), mc(1), 'filled', 'MarkerFaceColor', mc(2), 'MarkerEdgeColor','none');
                sc.MarkerFaceAlpha = modelEyeAlpha(pp);
            otherwise
                sc = scatter(imagePoints(idx,1), imagePoints(idx,2), mc(1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor',mc(2));
                sc.MarkerEdgeAlpha = modelEyeAlpha(pp);
        end
        % If we are plotting the pupil perimeter points, see if we would
        % like to label them
        if p.Results.showPupilTextLabels && strcmp(p.Results.modelEyeLabelNames{pp},'pupilPerimeterFront')
            % Put text labels for the pupil perimeter points so that we can follow
            % them through rotations and translations to validate the projection
            % model
            text(imagePoints(idx,1), imagePoints(idx,2), num2str(find(idx)));
        end

    end
end % loop over label names
hold off

% Get the rendered frame
renderedFrame=getframe(gcf);

end
