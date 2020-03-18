%% DEMO_ListingsLaw

% Create an eye with no corneal rotation, as this rotation can complicate
% the demo
sceneGeometry=createSceneGeometry('measuredCornealCurvature',[44.2410 45.6302 0 0 0]);

sceneGeometry.eye.stop.eccenFcnString = '@(x) 0';

modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

eyePoses=[-40 40 0 3; 0 40 0 3; 40 40 0 3; -40 0 0 3; 0 0 0 3; 40 0 0 3; -40 -40 0 3; 0 -40 0 3; 40 -40 0 3 ];
obey = [true,false];
figureNames = {'Obeys Listing''s Law','Disobeys Listing''s Law'};

for ff = 1:2
figure('Name',figureNames{ff});

for pp = 1:9
    subplot(3,3,pp)
    
    % Obtain the rendering of the model for this pose
    [figHandle, ~, renderedFrame] = renderEyePose(eyePoses(pp,:), sceneGeometry,...
        'visible', false, ...
        'nStopPerimPoints',8, ...
        'addPseudoTorsion',obey(ff),...
        'modelEyeAlpha',0.5,...
        'modelEyeLabelNames',modelEyeLabelNames,...
        'modelEyePlotColors',modelEyePlotColors);

    % Close the fig handle, as we will be displaying a mosaic of the
    % rendered images
    close(figHandle);
    
    % Display this eye render and clean up the plot
    subplot(3,3,pp);
    imshow(renderedFrame.cdata, 'Border', 'tight');
    axis off
    axis equal
    xlim([0 620]);
    ylim([0 480]);
    title(num2str(eyePoses(pose,:)));
    hold on
    
    % Super-impose a cross that shows the left-right / top-bottom
    % orientation of the pupil
    [~, ~, imagePoints, ~, ~, ~, pointLabels] = ...
        projectModelEye(eyePoses(pp,:), sceneGeometry,...
        'addPseudoTorsion',obey(ff),...
        'nStopPerimPoints',8);
    pupilPoints = find(strcmp(pointLabels,'pupilPerimeter'));
    imageSize = sceneGeometry.cameraIntrinsic.sensorResolution;
    A = imagePoints(pupilPoints(1),:);
    B = imagePoints(pupilPoints(5),:);
    plot([A(1) B(1)],[A(2),B(2)],'-r','LineWidth',3);
    A = imagePoints(pupilPoints(3),:);
    B = imagePoints(pupilPoints(7),:);
    plot([A(1) B(1)],[A(2),B(2)],'-r','LineWidth',3);
    
end
end
