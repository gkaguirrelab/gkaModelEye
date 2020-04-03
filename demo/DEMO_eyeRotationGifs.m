%% DEMO_eyeRotationGifs
% Create some gifs of the rotating eye
%
% Description:
%   Makes use of the Matlab Central function "gif"
%

% Save location
gifSaveRoot = '~/Desktop/EyeRotation';

% The default sceneGeometry
sg = createSceneGeometry;

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];

% File suffixes
fileSuffix = {'Azi','Ele','AziEle'};


%% Loop over gifs to create
for gg = 1:3
    
    % Construct the file name to save
    outputFile = [gifSaveRoot fileSuffix{gg} '.gif'];

    
    %% Loop over eyePoses
    for ii = 1:length(rotationValues)
        
        switch gg
            case 1
                eyePose = [rotationValues(ii) 0 0 3];
            case 2
                eyePose = [0 rotationValues(ii) 0 3];
            case 3
                eyePose = [rotationValues(ii) rotationValues(ii) 0 3];
        end
        
        if ii == 1
            [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', true, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
                'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});

            % This command opens the gif object           
            gif(outputFile);
            
        else
            delete(plotHandles(1:end))
            [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', false, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
                'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});
        end
        
        % This updates the gif
        gif

    end
end

% Close any open windows
close all
