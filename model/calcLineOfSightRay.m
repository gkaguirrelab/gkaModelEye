function [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
% Returns the path of the line of sight for a model eye
%
% Syntax:
%  [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
%
% Description
%   Given sceneGeometry, the routine identifies the ray that originates at
%   the fixation point, passes through the center of the entrance pupil,
%   and arrives at the foeca.
%
%   The radius of the aperture stop is assumed to be 3 mm unless set. The
%   fixation target is assumed to 1500 mm unless set.
%
%   The routine requires that the field sceneGeometry.eye.landmarks.fovea
%   be defined.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%   fixTargetDistance     - Scalar that is the Euclidean distance in mm of 
%                           the fixation target from the origin of the
%                           world coordinate frame. Defaults to 100.
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
%   fixationEyePose       - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           This is the pose of the eye for which the line
%                           of sight axis intersects the fixationTarget.
%   foveaDistanceError    - Scalar. The Euclidean distance in mm fromt the
%                           fovea to the point of the ray with the retina.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    [outputRayLoS,rayPathLoS,fixationEyePose]=calcLineOfSightRay(sceneGeometry,0.5);
    [outputRayVis,rayPathVis]=calcNodalRay(sceneGeometry.eye,sceneGeometry.eye.landmarks.fovea.geodetic);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPathVis,'outputRay',outputRayVis);
    plotOpticalSystem('newFigure',false,'rayPath',rayPathLoS,'outputRay',outputRayLoS,'rayColor','green');
    fprintf('Angle of the line-of-sight axis w.r.t. the optical axis:\n')
    fixationEyePose(1:2)
%}

% Parse inputs
if nargin==1
    stopRadius = 3;
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
fixTarget = @(theta,phi) [cosd(phi)*cosd(theta);cosd(phi)*sind(theta);sind(phi)].*fixTargetDistance;

% Anonymous function to return the Euclidean distance between the fovea and
% the point of intersection on the retina of the candidate line of sight
% ray
foveaDistance = @(coord) sqrt(sum((coord - sceneGeometry.eye.landmarks.fovea.coords).^2));

% Anonymous function to return the distance error of the line of sight
% intersecting the fixation target
myObj = @(p) foveaDistance(evalCandidateLineOfSight(sceneGeometry,stopRadius,fixTarget(p(1),p(2))));

% Define some options for the fmincon call
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% Define x0 and landmarks
x0 = sceneGeometry.eye.landmarks.fovea.degField(1:2);
lb = [-20 -20];
ub = [20 20];

% Search for the eyePose that results in the line of sight intersecting the
% fixation target
[p, foveaDistanceError] = fmincon(myObj,x0,[],[],[],[],lb,ub,[],opts);

% Get the angles with which the line of sight ray departs the fovea
[~,incidentRay] = evalCandidateLineOfSight(sceneGeometry,stopRadius,fixTarget(p(1),p(2)));
initialRay(:,1) =  incidentRay(:,1);
initialRay(:,2) = -incidentRay(:,2);

% Calculate and save the line of sight outputRay and the raypath
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;
[outputRay,rayPath] = rayTraceQuadrics(initialRay, opticalSystem);

% Calculate the fixationEyePose
[azimuth, elevation] = quadric.rayToAngles(outputRay);
fixationEyePose = [azimuth, elevation, 0, stopRadius];


end


%% LOCAL FUNCTIONS


%% evalCandidateLineOfSight
% Given a sceneGeometry, stopRadius, and a fixation coordinate in eyeWorld
% space, this routine traces a ray that originates at the fixation
% coordinate and strikes the center of the entrance pupil. The final point
% of intersection of the ray upon the retina is returned.
function [retinaCoords,outputRay] = evalCandidateLineOfSight(sceneGeometry,stopRadius,fixTarget)

% Place the camera at the fixation target
sceneGeometry.cameraPosition.translation = fixTarget([2 3 1]);

% Obtain the center of the entrance pupil for this eye pose
[~, ~, ~, ~, eyePoints, pointLabels] = pupilProjection_fwd([0 0 0 stopRadius], sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:));

% Find the ray that leaves the camera and strikes the center of the
% entrance pupil
R = quadric.normalizeRay([fixTarget';entrancePupilCenter-fixTarget']');

% Ray trace from the camera to the retina
outputRay = rayTraceQuadrics(R, sceneGeometry.refraction.cameraToRetina.opticalSystem);
retinaCoords = outputRay(:,1)';

end



