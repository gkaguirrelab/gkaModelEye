%% t02_pupilWithRotation
% Demonstrate that the pupil  moves differently from the eyeball due to the
% refractive effects of the cornea
%{
02: The refractive effect of the cornea varies with the relative position of the observer, causing the pupil to move differently from the sclera when the eye is in motion.
%}

% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = '~/Desktop/t02_pupilWithRotation.gif';

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];

%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [rotationValues(ii) -5 0 3];
    
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', true, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
        
        % This command opens the gif object
        gif(gifSaveName);
        
    else
        delete(plotHandles(1:end))
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
    end
    
    % This updates the gif
    gif
    
end

% Close any open windows
close all
