%% DEMO_ListingsLaw
% Demonstrate the application of psuedo torsion to the eyePose
%
% Description
%   Head-fixed eye movements generally obey Listing's Law, which observes
%   that there is no change in the "true" torsion of the eye (with respect
%   to its optical axis) with a change in position. The current model
%   describes eye rotations using "Fick coordinates". To create eye
%   movements that obey Listing's Law, it is necessary to
%   add "pseudo torsion" to the eyePose. This correction takes place
%   intenally within the the projectModelEye.m function, with a call to the
%   sub-routine addPseudoTorsion.m.
%
%   The addition of pseudo torsion is under the control of the
%   addPseudoTorsion key, which is set to true by default. This demo
%   illustrates the appearance of the eye in different poses with and
%   without the addition of pseudo torsion.
%
%   The calculation of pseudoTorsion is performed relative to the "primary
%   position" of the eye, which is the position from which any other eye
%   pose can be achieved by rotation of the eye about a single axis. This
%   value is set in:
%
%       sceneGeometry.eye.rotationCenters.primaryPosition
%
%   The default value for the primary position is [0,0], and this value is
%   defined within the camera coordinate space, in which [0, 0] indicates
%   that the optical axis of the eye is aligned with the optical axis of
%   the camera. There is no requirement in practice that the fixation point
%   of the eye, the primary position of the eye, and the angle of these
%   points with respect to the optical axis of the camera all align.
%


%% Create the sceneGeometry

% Create a sceneGeometry in which the cornea is aligned with the optical
% axis of the eye. Forcing this alignment makes the demo a bit clearer.
sceneGeometry=createSceneGeometry('measuredCornealCurvature',[44.2410 45.6302 0 0 0]);

% For the same reason, force the aperture stop of the iris to be perfectly
% circular
sceneGeometry.eye.stop.eccenFcnString = '@(x) 0';


%% Set up variables for the plots

modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};
eyePoses=[-40 40 0 3; 0 40 0 3; 40 40 0 3; -40 0 0 3; 0 0 0 3; 40 0 0 3; -40 -40 0 3; 0 -40 0 3; 40 -40 0 3 ];
obey = [true,false];
figureNames = {'Obeys Listing''s Law','Disobeys Listing''s Law'};


%% Make two figures
for ff = 1:2
    
    % Create the figure
    figure('Name',figureNames{ff});
    
    % Loop over eye poses
    for pp = 1:size(eyePoses,1)
        
        % Into the subplot you go
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
        title(num2str(eyePoses(pp,:)));
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
