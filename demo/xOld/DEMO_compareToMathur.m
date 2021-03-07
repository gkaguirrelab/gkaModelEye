%% DEMO_compareToMathur
% Compare the shape of the modeled pupil ellipse to empirical measurements
%
% Description:
%   An empirical measurement of the aspect ratio of the pupil was made
%   with the camera positioned at a range of angles with respect to the
%   line of sight of the eye:
%
%       Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil shape
%       as viewed along the horizontal visual field." Journal of vision
%       13.6 (2013): 3-3.
%
%   This demo compares the output of the modeled pupil with the measured
%   shape.


% The refractive error of the subject for the average Mathur data.
sphericalAmetropia = (6*1.2-11*2.9)/30;

% The Mathur subjects were  dilated with 1% Cyclopentolate which in adults
% produces an entrance pupil of ~6 mm:
%
%   Kyei, Samuel, et al. "Onset and duration of cycloplegic action of 1%
%   Cyclopentolate-1% Tropicamide combination." African health sciences
%   17.3 (2017): 923-932.
%
entrancePupilDiam = 6;

% Calculate the aperture stop diameter that corresponds to this
% entrance pupil diameter.
%{
    entranceRadius = 6/2;
    % Prepare scene geometry and eye pose aligned with visual axis
    sceneGeometry = createSceneGeometry();
    % Obtain the pupil area in the image for the entrance radius
    % assuming no ray tracing
    sceneGeometry.refraction = [];
    pupilImage = projectModelEye([0, 0, 0, entranceRadius],sceneGeometry);
    stopArea = pupilImage(3);
    % Add the ray tracing function to the sceneGeometry
    sceneGeometry = createSceneGeometry();
    % Search across stop radii to find the value that matches the observed
    % entrance area.
    myPupilEllipse = @(radius) projectModelEye([0, 0, 0, radius],sceneGeometry);
    myArea = @(ellipseParams) ellipseParams(3);
    myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
    stopRadius = fminunc(myObj, entranceRadius)
%}
stopDiam = 2.6547*2;


% Subjects in the Mathur study fixated a point 3 meters away
targetDistance = 3000;
accommodation = 1000/3000;

% Create an eye
eye = modelEyeParameters('sphericalAmetropia',sphericalAmetropia,...
    'accommodation',accommodation,...
    'spectralDomain','vis');

% Obtain the eyePose for fixating the target (which is at the location of
% camera position [0,0].
fixationEyePose=calcFixationPose(eye,[0 0],targetDistance,true,stopDiam/2);

% Set the rotation center of the eye to the origin of the longitudinal
% axis, which is also the apex of the un-rotated cornea. This allows us to
% simulate rotating a camera around this point.
eye.rotationCenters.azi = [0 0 0];
eye.rotationCenters.ele = [0 0 0];

% Put that eye in a sceneGeometry
sceneGeometry = createSceneGeometry(...
    'eye',eye,...
    'cameraTranslation',[0; 0; 100],...
    'spectralDomain','vis');

% The range for our plots
viewingAngleDeg = -80:5:70;

% Loop over the viewing angles and calculate the diameter ratio
for vv = 1:length(viewingAngleDeg)    
    
    % Our model rotates the eye. For the right eye, a positive azimuth rotates
    % the eye such that the center of the pupil moves to the right of the
    % image. This means that a positive azimuth corresponds to the camera being
    % positioned in the temporal portion of the visual field. So, we must sign
    % reverse the interpretation of our azimuthal values for measurements made
    % in the right eye to correspond to the Mathur results. Additionally, we
    % need to adjust for alpha: the angle between the pupil and visual axes of
    % the eye. The coordinates of our model eye are based around the pupil
    % axis. Therfore, we need to calculate a rotation that accounts for the
    % Mathur viewing angle and alpha.
    azimuthDeg = (-viewingAngleDeg(vv))+fixationEyePose(1);
    elevationDeg = zeros(size(viewingAngleDeg(vv)))+fixationEyePose(2);
    
    % Assemble the eyePose
    eyePose=[azimuthDeg elevationDeg 0 stopDiam/2];
    
    % First, perform the forward projection to determine where the center of
    % the entrance pupil is located in the sceneWorld coordinates
    % Obtain the center of the entrance pupil for this eye pose.
    [~, ~, ~, worldPoints, ~, ~, pointLabels] = projectModelEye(eyePose, sceneGeometry, 'nStopPerimPoints', 16);
    pupilCenter = nanmean(worldPoints(strcmp(pointLabels,'pupilPerimeter'),:));
    
    % We will translate the camera to be centered on the entrance pupil. This
    % is an attempt to match the arrangement of the Mathur study, in which the
    % examiner adjusted the camera to be centered on the pupil.
    cameraTrans = pupilCenter';
    
    % Now, measure the pupil diameter ratio
    pupilEllipseOnImagePlane = ...
        projectModelEye(eyePose, sceneGeometry,'cameraTrans',cameraTrans,...
        'nStopPerimPoints',16,'replaceReflectedPoints',false);
    
    % Obtain the ellipse in explicit format
    p = ellipse_transparent2ex(pupilEllipseOnImagePlane);
    
    % Calculate the diameter ratio, and account for the possibility of an
    % ellipse that has a horizontal major axis
    theta = pupilEllipseOnImagePlane(5);
    if theta < pi/4
        diamRatio(vv)=p(3)./p(4);
    else
        diamRatio(vv)=p(4)./p(3);
    end
    
end

% This is Eq 9 from Mathur 2013, which specifies the horizontal to vertical
% ratio of the entrance pupil from different viewing angles relative to
% fixation
mathurEq9 = @(viewingAngleDeg) 0.99.*cosd((viewingAngleDeg+5.3)/1.121);

% Plot the results.
figHandle1 = figure();

plot(viewingAngleDeg,mathurEq9(viewingAngleDeg) ,'-','Color',[.5 .5 .5]);
hold on
plot(viewingAngleDeg,diamRatio,'-','Color',[1 0 0]);

pbaspect([1 1.5 1])
xlim([-90 90]);
xticks([-75 -50 -25 0 25 50 75])
ylim([0 1]);
xlabel('Viewing angle [deg]')
ylabel('Pupil Diameter Ratio')
legend({'Mathur 2013','Current model'},'Location','southeast');


