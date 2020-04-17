%% t06_fixationPose
% Demonstrate the appearance of eyes that are fixated upon a distant point
%{
06: As the fovea is typically located inferotemporally on the retina, when fixating a distant target, the eyes are turned slightly down and out.
%}

% Right eye, left eye
laterality = {'right','left'};

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% Loop over eyes
for ii = 1:length(laterality)
    
    % Create the sceneGeometry
    sceneGeometry=createSceneGeometry('calcLandmarkFovea',true,'eyeLaterality',laterality{ii});
    
    % Define a stop radius for the eye. This value produces a pupil that is
    % about 3.5 mm in diamter
    stopRadius = 1.53;
    
    % Set the target distance, which is on the optical axis of each eye. We
    % want this far enough away that we are not modeling vergence of the
    % eyes.
    fixTargetDistance = 1500;
    
    % Get the eyePose that places the fixation target on the line of sight
    [~,~,fixEyePose] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance);
    
    % Render the eye with corneal refraction
    [~, ~, frame] = renderEyePose(fixEyePose, sceneGeometry, ...
        'visible', false, ...
        'modelEyeLabelNames',modelEyeLabelNames,...
        'modelEyePlotColors',modelEyePlotColors);
    
    % Save this frame
    framesToMontage(:,:,:,ii) = frame.cdata;

end

% Set up a figure
figure

% Show the montage
montage(framesToMontage);

