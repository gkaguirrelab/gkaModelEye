%% t07_torsion
% Demonstrate the appearance of an eye undergoing torsion
%{
07: In addition to horizontal and vertical rotation, the eye can undergo rotation about the optical axis, which is termed "torsion". What is the torsion of the eye when it is gazing at a particular point in space? 
%}

% Create the sceneGeometry
sceneGeometry=createSceneGeometry();

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = '~/Desktop/t07_torsion.gif';

% The angles across which the eye will rotate
rotationValues = [0:1:10,10:-1:-10,-10:1:-1];


%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [0 0 rotationValues(ii) 3];
    
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
