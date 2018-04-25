function [figHandle, renderedFrame] = renderEyePose(eyePose, sceneGeometry, varargin)
% Creates an image of the eye for a given eyePose
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
%  'newFigure'            - Logical. Determines if we create a new figure.
%
% Outputs:
%   figHandle             - Handle to a created figure. Empty if a new
%                           figure was not requested.
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
    figure
    for azi = -35:35:35
        for ele = -35:35:35
            eyePose = [azi ele 0 1];
            renderEyePose(eyePose,sceneGeometry,'newFigure',false,'modelEyeLabelNames',{'pupilPerimeter'},'modelEyePlotColors',{'.g'});
        end
    end
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eyePose',@(x)(isnumeric(x) && all(size(x)==[1 4])));
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('newFigure',true,@islogical);
p.addParameter('modelEyeLabelNames', {'aziRotationCenter', 'eleRotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex'}, @iscell);
p.addParameter('modelEyePlotColors', {'>r' '^m' '.w' '.b' '*g' '.y' '*y'}, @iscell);


% parse
p.parse(eyePose, sceneGeometry, varargin{:})

% Grab the image size
imageSizeX = sceneGeometry.cameraIntrinsic.sensorResolution(1);
imageSizeY = sceneGeometry.cameraIntrinsic.sensorResolution(2);

% A blank frame to initialize each frame
blankFrame = zeros(imageSizeY,imageSizeX)+0.5;

% Open a figure
if p.Results.newFigure
    figHandle = figure;
    imshow(blankFrame, 'Border', 'tight');
end

% Prepare the figure
hold on
axis off
axis equal
xlim([0 imageSizeX]);
ylim([0 imageSizeY]);

% Obtain the pupilProjection of the model eye to the image plane
[pupilEllipseParams, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd(eyePose, sceneGeometry, 'fullEyeModelFlag', true, 'nIrisPerimPoints',20, 'removeOccultedPoints', false);

% Loop through the point labels present in the eye model
for pp = 1:length(p.Results.modelEyeLabelNames)
    idx = strcmp(pointLabels,p.Results.modelEyeLabelNames{pp});
    if strcmp(p.Results.modelEyeLabelNames{pp},'pupilPerimeter')
        % Just before we plot the pupil perimeter points, add the
        % pupil fit ellipse
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilEllipseParams));
        fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
        % superimpose the ellipse using fimplicit or ezplot (ezplot
        % is the fallback option for older Matlab versions)
        if exist('fimplicit','file')==2
            fimplicit(fh,[1, imageSizeX, 1, imageSizeY],'Color', 'g','LineWidth',1);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis off;
        else
            plotHandle=ezplot(fh,[1, imageSizeX, 1, imageSizeY]);
            set(plotHandle, 'Color', p.Results.pupilColor)
            set(plotHandle,'LineWidth',1);
        end
    end
    plot(imagePoints(idx,1), imagePoints(idx,2), p.Results.modelEyePlotColors{pp})
end
hold off

% Get the rendered frame
renderedFrame=getframe(gcf);

end
