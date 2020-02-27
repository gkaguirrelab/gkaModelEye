%% DEMO_correctiveLenses
% Demonstrate the effect of corrective lenses in the model
%
% Description:
%   Data is sometimes collected from subjects while they are wearing either
%   contact lenses or spectacles. This routine demonstrates the change in
%   the appearance of the pupil in the image plane as a consequence of the
%   presence of corrective lenses.
%

clear all
close all
clc


% Define some variables for plotting the model eye
diopterColors = {'g','r','b'};
blankFrame = zeros(480,640)+0.5;


%% Present Figure 1 -- Spectacle lens virtual image displacement

fprintf(['Figure 1 shows the outline of the pupil on the image plane for an\n' ...
    'emmetropic eye in three different azimuthal poses. The pupil ellipse is drawn in \n' ...
    'green for a projection that does not include a corrective lens, red for \n' ...
    'the projection with a -4 diopter spectacle lens, and blue for a +4 diopter \n' ...
    'lens. As can be seen, a spectacle lens shifts the virtual image of the \n'...
    'pupil. The corresponding glint locations are also indicated. \n' ...
    'The left/right asymmetry for + and - azimuth is due to the tilt \n' ...
    'of the corneal ellipse. \n\n']);

% Open a figure and define the eye poses and refractive corrections
figure(1)
eyePoses=[-20 0 0 3; 0 0 0 3; 20 0 0 3];
lensRefractionDiopters = [0, -4, 4];

for dd = 1:length(lensRefractionDiopters)
    
    % Obtain the sceneGeometry and ray tracing functions
    sceneGeometry = createSceneGeometry('sphericalAmetropia',0,'spectacleLens',lensRefractionDiopters(dd));
    sceneGeometry.cameraPosition.translation(3)=60;
    
    for pose = 1:size(eyePoses,1)
        % Perform the projection and request the full eye model
        [pupilEllipseOnImagePlane, glintCoord, imagePoints, ~, ~, ~, pointLabels] = projectModelEye(eyePoses(pose,:),sceneGeometry,'nStopPerimPoints',8);
        % plot the pupil ellipse
        subplot(1,3,pose);
        if dd==1
            imshow(blankFrame, 'Border', 'tight');
        hold on
        axis off
        axis equal
        xlim([0 640]);
        ylim([0 480]);
        end
        % Plot the pupil ellipse
        addTransparentEllipseToFigure(pupilEllipseOnImagePlane,640,480,diopterColors{dd});
        plot(glintCoord(1),glintCoord(2),['*' diopterColors{dd}]);
        title(num2str(eyePoses(pose,:)));
    end
    drawnow
end

