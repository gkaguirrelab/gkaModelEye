%% t02_pupilWithRotation
% Demonstrate that the pupil  moves differently from the eyeball due to the
% refractive effects of the cornea


% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = '~/Desktop/EyeRotation.gif';

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];


%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [rotationValues(ii) -5 0 3];
    
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', true, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y'});
        
        % This command opens the gif object
        gif(gifSaveName);
        
    else
        delete(plotHandles(1:end))
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y'});
    end
    
    % This updates the gif
    gif
    
end

% Close any open windows
close all
