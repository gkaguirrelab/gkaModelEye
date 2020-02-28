function [outputRay,rayPath,fixEyePose,fixTargetCoords,foveaDistanceError] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
% Returns the path of the line of sight for a model eye
%
% Syntax:
%  [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
%
% Description
%   Given sceneGeometry, the routine identifies the ray that originates at
%   the fixation point, passes through the center of the entrance pupil,
%   and arrives at the fovea.
%
%   This ray is defined as the "line of sight" for the eye.
%
%   If not defined, the radius of the aperture stop is set to provide an
%   entrance pupil diameter of ~2 mm, which empirically produces the
%   highest degree of acuity in normal observers. The fixation target is
%   assumed to 1500 mm unless set.
%
%   The routine requires that the field sceneGeometry.eye.landmarks.fovea
%   be defined.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   stopRadius            - Scalar. The radius of the aperture stop
%   fixTargetDistance     - Scalar that is the Euclidean distance in mm of 
%                           the fixation target from the origin of the
%                           world coordinate frame. Defaults to 1500.
%
% Outputs:
%   outputRay             - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   fixEyePose            - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           This is the pose of the eye for which the line
%                           of sight axis intersects the fixationTarget.
%   fixTargetCoords       - A 3x1 vector that gives the location of the
%                           fixation target in the world coordinate space
%                           in units of mm.
%   foveaDistanceError    - Scalar. The Euclidean distance in mm fromt the
%                           fovea to the point of the ray with the retina.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    [outputRayLoS,rayPathLoS,fixationEyePose]=calcLineOfSightRay(sceneGeometry);
    [outputRayVis,rayPathVis]=calcNodalRay(sceneGeometry.eye,sceneGeometry.eye.landmarks.fovea.geodetic);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPathVis,'outputRay',outputRayVis);
    plotOpticalSystem('newFigure',false,'rayPath',rayPathLoS,'outputRay',outputRayLoS,'rayColor','green');
    outline = sprintf('Angle of the line-of-sight axis w.r.t. the optical axis: [%2.2f, %2.2f]\n',fixationEyePose(1:2));
    fprintf(outline);
%}
%{
    % Derive the axial length of the eye along the line-of-sight axis
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    [~,rayPathLoS]=calcLineOfSightRay(sceneGeometry);
    al = sqrt(sum((rayPathLoS(:,1)-rayPathLoS(:,end)).^2));
    outline = sprintf('The axial length of the eye along the line-of-sight axis is: %2.2f mm\n',al);
    fprintf(outline);
%}

% Code to determine the stop radius that corresponds to a pupil diameter of
% 2 mm. This value is used as it is found to provide peak acuity for normal
% observers.
%{
    entranceRadius = 2/2;
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
    stopRadius = fminunc(myObj, entranceRadius);
    outline = sprintf('A 2mm entrance pupil corresponds to a %2.2fmm stop radius\n',stopRadius);
    fprintf(outline);
%}
    
% Parse inputs
if nargin==1
    stopRadius = 0.8693;
    fixTargetDistance = 1500;
end

if nargin==2
    fixTargetDistance = 1500;
end

% Check that the sceneGeometry eye has a foveal landmark
if ~isfield(sceneGeometry.eye.landmarks,'fovea')
    error('The sceneGeometry does not have a fovea defined. Set calcLandmarkFovea to true in the call to createSceneGeometry')
end

% Anonymous function to return a fixation coordinate at a given azimuth and
% elevation relative to the origin of the eye coordinate space (the
% un-rotated corneal apex).
% -180 <= theta <= 180 (azimuth)
% -90 <= phi <= 90 (elevation)
fixEyeWorldTargetFunc = @(theta,phi) [cosd(phi)*cosd(theta);cosd(phi)*sind(theta);sind(phi)].*fixTargetDistance;

% Anonymous function to return the Euclidean distance between the fovea and
% the point of intersection on the retina of the candidate line of sight
% ray
foveaDistance = @(coord) sqrt(sum((coord - sceneGeometry.eye.landmarks.fovea.coords).^2));

% Anonymous function to return the distance error of the line of sight
% intersecting the fixation target
myObj = @(p) foveaDistance(evalCandidateLineOfSight(sceneGeometry,stopRadius,fixEyeWorldTargetFunc(p(1),p(2))));

% Define x0 and landmarks
x0 = sceneGeometry.eye.landmarks.fovea.degField(1:2);

% define some search options
options = optimset('Display','off');

% Search for the eyePose that results in the line of sight intersecting the
% fixation target
[p, foveaDistanceError] = fminsearch(myObj,x0,options);

% Obtain and save the fixation target coords
fixTargetEyeWorldCoords = fixEyeWorldTargetFunc(p(1),p(2));

% Get the angles with which the line of sight ray departs the fovea
[~,incidentRay] = evalCandidateLineOfSight(sceneGeometry,stopRadius,fixTargetEyeWorldCoords);
initialRay(:,1) =  incidentRay(:,1);
initialRay(:,2) = -incidentRay(:,2);

% Calculate and save the line of sight outputRay and the raypath
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;
[outputRay,rayPath] = rayTraceQuadrics(initialRay, opticalSystem);

% Calculate the fixationEyePose
[azimuth, elevation] = quadric.rayToAngles(outputRay);
fixEyePose = [azimuth, elevation, 0, stopRadius];

% Re-arrrange the eyeWorld coordinates to return the target in World coords
fixTargetCoords = fixTargetEyeWorldCoords([2 3 1]);


end


%% LOCAL FUNCTIONS


%% evalCandidateLineOfSight
% Given a sceneGeometry, stopRadius, and a fixation coordinate in eyeWorld
% space, this routine traces a ray that originates at the fixation
% coordinate and strikes the center of the entrance pupil. The final point
% of intersection of the ray upon the retina is returned.
function [retinaCoords,outputRay] = evalCandidateLineOfSight(sceneGeometry,stopRadius,fixEyeWorldTarget)

% Place the camera at the fixation target
sceneGeometry.cameraPosition.translation = fixEyeWorldTarget([2 3 1]);

% Obtain the center of the entrance pupil for this eye pose
[~, ~, ~, ~, ~, eyePoints, pointLabels] = projectModelEye([0 0 0 stopRadius], sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:));

% Find the ray that leaves the camera and strikes the center of the
% entrance pupil
R = quadric.normalizeRay([fixEyeWorldTarget';entrancePupilCenter-fixEyeWorldTarget']');

% Ray trace from the camera to the retina
outputRay = rayTraceQuadrics(R, sceneGeometry.refraction.cameraToRetina.opticalSystem);
retinaCoords = outputRay(:,1)';

end



