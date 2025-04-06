%% t08_ListingsLaw
% Demonstrate the appearance of a moving eye that obeys Listing's Law
%{
08: For voluntary, head-fixed eye movements, Listing's Law specifies eye torsion for any pose by reference to its "primary position". The model on the right uses "Fick" eye rotation coordinates without "pseudo torsion" correction.
%}
% An understanding of this topic requires grappling with different
% coordinate systems for specifying 3D rotations. For details see:
%   project/addPseudoTorsion.m
%

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = 'demo/t08_ListingsLaw.gif';

% Create an eye with a specific measuredCornealCurvature sceneGeometry in which the cornea is aligned with the optical
% axis of the eye. Forcing this alignment makes the demo a bit clearer.
sceneGeometry=createSceneGeometry('kvals',[44.2410 45.6302 0 0 0]);

% For the same reason, force the aperture stop of the iris to be perfectly
% circular
sceneGeometry.eye.stop.eccenFcnString = '@(x) 0';

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y'};

% The angles across which the eye will rotate
rotationValues = [0:1:30,30:-1:-30,-30:1:-1];

% Plot labels
plotLabels = {'Obeys Listing''s Law','Disobeys Listing''s Law'};

% Obedience to Listing's Law
obedience = [true,false];

% Create a figure
figHandle = figure();


%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    % Set the eyePose for this frame
    eyePose = [rotationValues(ii) rotationValues(ii) 0 3];
        
    %% Render following and violating Listing's Law
    % We tell the model to "addPseudoTorsion" or not so that they eye pose
    % is in agreement with Listing's Law.
    
    for ff = 1:2
        
        % Obtain the rendering of the model for this pose
        tmpFigHandle = renderEyePose(eyePose, sceneGeometry,...
            'visible', false, ...
            'nStopPerimPoints',8, ...
            'addPseudoTorsion',obedience(ff),...
            'modelEyeAlpha',0.5,...
            'modelEyeLabelNames',modelEyeLabelNames,...
            'modelEyePlotColors',modelEyePlotColors);
        
        % Super-impose a cross that shows the left-right / top-bottom
        % orientation of the pupil
        [~, ~, imagePoints, ~, ~, ~, pointLabels] = ...
            projectModelEye(eyePose, sceneGeometry,...
            'addPseudoTorsion',obedience(ff),...
            'nStopPerimPoints',8);
        pupilPoints = find(strcmp(pointLabels,'pupilPerimeter'));
        A = imagePoints(pupilPoints(1),:);
        B = imagePoints(pupilPoints(5),:);

        set(0, 'CurrentFigure', tmpFigHandle)
        hold on
        plot([A(1) B(1)],[A(2),B(2)],'-r','LineWidth',3);
        A = imagePoints(pupilPoints(3),:);
        B = imagePoints(pupilPoints(7),:);
        plot([A(1) B(1)],[A(2),B(2)],'-r','LineWidth',3);
        
        % Get the frame
        drawnow;
        frame=getframe(tmpFigHandle);
        
        % Close the invisible figure
        close(tmpFigHandle);
        
        % Add a text label to the frame
        frame.cdata = insertText(frame.cdata,[20 20],plotLabels{ff},'FontSize',20);
        
        % Store the frame
        framesToMontage(:,:,:,ff) = frame.cdata;
        
    end
    
    
    %% Display the montage
    
    % Create the montage
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

close(figHandle)

