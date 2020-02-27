

sg = createSceneGeometry;

% These are the median (across subject) rotation centers from the TOME
% study
sg.eye.rotationCenters.azi(1:2) = [-12.15, 0.64];
sg.eye.rotationCenters.ele(1:2) = [-10.82, 0.00];

rotationValues = [0:1:15,15:-1:-15,-15:1:-1];

% Create and save the example of azimuthal rotation
outputFile = '~/Desktop/eyeRotationAzi.gif';
for ii = 1:length(rotationValues)
    eyePose = [rotationValues(ii) 0 0 3];
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', true, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});
        gif(outputFile);
    else
        delete(plotHandles(2:end))
        [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', false, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});
    end
    gif
end


% Create and save the example of elevational rotation
outputFile = '~/Desktop/eyeRotationEle.gif';
for ii = 1:length(rotationValues)
    eyePose = [0 rotationValues(ii) 0 3];
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', true, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});
        gif(outputFile);
    else
        delete(plotHandles(2:end))
        [~, plotHandles] = renderEyePose(eyePose, sg, 'newFigure', false, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01'}, ...
            'modelEyePlotColors', {'.w' 'ob' '-g' '.y' '*r'});
    end
    gif
end

close all
