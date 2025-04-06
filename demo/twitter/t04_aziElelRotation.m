%% t04_aziEleRotation
% Demonstrate an exagerated version of the appearance of the eye for
% horizontal and vertical rotation. This routine creates a gif file.
%{
04: The effect of different effective centers for horizontal and vertical rotation is quite small and usually ignored, but is here exaggerated by 1.2x in the vertical, illustrating a "tilt" component of vertical motion.
%}

% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% Exagerate the shallow position of the elevational rotation center by 1.2x
sceneGeometry.eye.rotationCenters.ele = sceneGeometry.eye.rotationCenters.ele ./ 1.2;

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = 'demo/t04_aziEleRotation.gif';

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];

% Model components to render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% Create a figure
figHandle = figure();

%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [rotationValues(ii) 0 0 3];
    [~, ~, frameAzi] = renderEyePose(eyePose, sceneGeometry, 'visible', false, ...
        'modelEyeLabelNames', modelEyeLabelNames , ...
        'modelEyePlotColors', modelEyePlotColors );
    
    eyePose = [0 rotationValues(ii) 0 3];
    [~, ~, frameEle] = renderEyePose(eyePose, sceneGeometry, 'visible', false, ...
        'modelEyeLabelNames', modelEyeLabelNames , ...
        'modelEyePlotColors', modelEyePlotColors );
    
    % Create the montage
    framesToMontage(:,:,:,1) = frameAzi.cdata;
    framesToMontage(:,:,:,2) = frameEle.cdata;
    figure(figHandle)
    montage(framesToMontage);
    drawnow
    
    if ii == 1
        % This command opens the gif object
        gif(gifSaveName,'frame',figHandle);
    else
        % This updates the gif
        gif
    end
    
    
end

% Close any open windows
close all
