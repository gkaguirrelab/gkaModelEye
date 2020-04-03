%% DEMO_eyePoseParams
% Demonstrate the interpretation of the eyePose values
%
% Description:
%   The routines express the pose of the eye in a 1x4 vector of eyePoses
%   of the form [azimuth, elevation, torsion, pupilRadius]
%
%   Here we demonstrate the interpretation of these parameters.
%


%% Housekeeping
clear vars
close all
clc


% Some info for the console
fprintf(['The pose of the eye is described by the parameters:\n\n' ...
    '\t[azimuth, elevation, torsion, pupilRadius]\n\n' ...
    'The three rotation variables are in units of degrees, and are in the\n' ...
    '"Fick" coordinates.\n' ...
    '(i.e., rotating) coordinate frame.\n\n']);



%% Present Figure 1

figure(1)
eyePoses=[-20 20 0 3; 0 20 0 3; 20 20 0 3; -20 0 0 3; 0 0 0 3; 20 0 0 3; -20 -20 0 3; 0 -20 0 3; 20 -20 0 3 ];

% Obtain the sceneGeometry
sceneGeometry = createSceneGeometry();

for pose = 1:size(eyePoses,1)
    % Obtain the rendering of the model for this pose
    [figHandle, ~, renderedFrame] = renderEyePose(eyePoses(pose,:), sceneGeometry, 'visible', false);
    % Close the fig handle, as we will be displaying a mosaic of the
    % rendered images
    close(figHandle);
    % plot
    subplot(3,3,pose);
    imshow(renderedFrame.cdata, 'Border', 'tight');
    axis off
    axis equal
    xlim([0 620]);
    ylim([0 480]);
    title(num2str(eyePoses(pose,:)));
end
fprintf(['Figure 1 shows the pose of the eye across positive and negative values\n' ...
    'of azimuth and elevation. While,the specified torsion of the eye is zero,\n' ...
    'the model eye obeys Listing''s Law, so the projection routine adds \n' ...
    '''pseudo torsion'' internally in the model to keep the actual torsion zero.\n\n']);


%% Present Figure 2
figure(2)
eyeSides = {'right','left'};
modelEyeLabelNames = {'irisPerimeter' 'pupilEllipse' 'pupilCenter' 'cornealApex'};
modelEyePlotColors = {'.b' '-g' '+g' '+y'};
partsToPlot = [3 4 7 8 9];


for laterality = 1:2
    % prepare the model eye for this laterality
    sceneGeometry = createSceneGeometry('eyeLaterality',eyeSides{laterality},'calcLandmarkFovea',true);
    sceneGeometry.cameraPosition.translation(3)=50;
    [figHandle, ~, renderedFrame] = renderEyePose([0 0 0 3], sceneGeometry, 'visible', false,'modelEyeLabelNames',modelEyeLabelNames,'modelEyePlotColors',modelEyePlotColors);
    % Close the fig handle, as we will be displaying a mosaic of the
    % rendered images
    close(figHandle);
    % plot
    subplot(1,2,laterality);
    imshow(renderedFrame.cdata, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 620]);
    ylim([0 480]);
    text(320,50,[eyeSides{laterality} ' eye'],'HorizontalAlignment','center');
    
    % Add a point that corresponds to the visual axis at the pupil
    % plane
    [~, ~, imagePointsFixationAxis, ~, ~, pointLabelsFixationAxis] = ...
        projectModelEye([sceneGeometry.eye.landmarks.fovea.degField(1) sceneGeometry.eye.landmarks.fovea.degField(2) 0 3],sceneGeometry,'fullEyeModelFlag',true);
    
    idx = strcmp(pointLabelsFixationAxis,'pupilCenter');
    plot(imagePointsFixationAxis(idx,1), imagePointsFixationAxis(idx,2), '+c')    
    hold off
end
fprintf(['Figure 2 top shows the perimeter of the pupil (green) and\n' ...
    'iris (blue) for eyePoses [0 0 0 3] for the right and left eye. The\n' ...
    'axis of the camera is aligned with the optical axis of the model eye. \n' ...
    'Note that the pupil is slight ellipitical, with a vertical axis. \n' ...
    'Additionally, the cornea is aspheric and rotated out of alignment \n'...
    'with the optical axis, which adds additional non-circularity to the pupil.\n' ...
    'The yellow plus is the corneal apex, which is displaced nasally\n\n']);

