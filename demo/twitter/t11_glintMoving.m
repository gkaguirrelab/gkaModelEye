%% t11_glintMoving
% Illustrate the ray path from light to camera to produce a glint
%{
11: The glint (*) shifts with eye movement, influenced by the non-sphericity of the cornea and the center of rotation of the eye. The position of the glint relative to the pupil center is often used in eye tracking to derive gaze position.
%}
% For details see:
%   Merchant, John, Richard Morrissette, and James L. Porterfield. "Remote
%   measurement of eye direction allowing subject motion over one cubic
%   foot of space." IEEE transactions on biomedical engineering 4 (1974):
%   309-317.
%


% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = '~/Desktop/t11_glintMoving.gif';

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y' '*r'};

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];

% Get the image dimensions
imDims = sceneGeometry.cameraIntrinsic.sensorResolution;

%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [rotationValues(ii) -5 0 2];
    
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', true, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
        
        % Add cross-hairs to assist in detection of the glint motion
        hold on
        plot([0 imDims(1)],[imDims(2)/2 imDims(2)/2],'-','Color',[0.75 0.75 0.75]);
        plot([imDims(1)/2 imDims(1)/2],[0 imDims(2)],'-','Color',[0.75 0.75 0.75]);
        
        % This command opens the gif object
        gif(gifSaveName);
        
    else
        delete(plotHandles(1:end))
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
    end
    
    % Add cross-hairs to assist in detection of the glint motion
    hold on
    plot([0 imDims(1)],[imDims(2)/2 imDims(2)/2],'-','Color',[0.75 0.75 0.75]);
    plot([imDims(1)/2 imDims(1)/2],[0 imDims(2)],'-','Color',[0.75 0.75 0.75]);
    
    % This updates the gif
    gif
    
end

% Close any open windows
close all

