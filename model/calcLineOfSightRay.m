function [outputRay,rayPath,distanceErrorEntrancePupil,distanceErrorFixationTarget] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
% Returns the path of the line of sight for a model eye
%
% Syntax:
%  [outputRay,rayPath] = calcLineOfSightRay(eye,G,X,cameraMedium)
%
% Description
%   Given sceneGeometry, the routine identifies the pose of the eye and a
%   ray, such that the ray originates at the fixation target, passes
%   through the center of the entrance pupil, and intersects the fovea.
%
%   The radius of the aperture stop is assumed to be 2mm unless set. The
%   fixation target is assumed to be the camera location unless set.
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
%   fixationPose          - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           This is the pose of the eye for which the line
%                           of sight axis intersects the fixationTarget.
%   distanceErrorEntrancePupil,distanceErrorFixationTarget - Scalars. The
%                           distance in mm that the ray passes by each of
%                           these targets.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    [outputRay,rayPath]=calcLineOfSightRay(sceneGeometry);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
    [angle_p1p2,angle_p1p3] = quadric.rayToAngles(outputRay);
%}

% Parse inputs
if nargin==1
    stopRadius = 2;
    fixTargetDistance = 1500;
end

if nargin==2
    fixTargetDistance = 1500;
end

% Check that the sceneGeometry eye has a foveal landmark
if ~isfield(sceneGeometry.eye.landmarks,'fovea')
    error('The sceneGeometry does not have a fovea defined. Call createSceneGeometry with calcLandmarkFovea set to true')
end

% Define the fixation target
fixTargetWorld = [0; 0; fixTargetDistance];

% Anonymous function to return the distance error of the line of sight
% intersecting the fixation target
myError = @(eyePose) evalCandidateLineOfSight(eyePose,sceneGeometry,fixTargetWorld);

% Define some options for the fmincon call
options = optimoptions(@fmincon,...
    'Display','off');

% Set the bounds on the eyePose search
x0 = [0 0 0 stopRadius];
lb = [-20 -10 0 stopRadius];
ub = [20 10 0 stopRadius];

% Search for the eyePose that results in the line of sight intersecting the
% fixation target
[fixationPose, distanceErrorFixationTarget] = fmincon(myError,x0,[],[],[],[],lb,ub,[],options);

% Obtain the initial angles with respect to the optical axis which the ray
% departed the fovea
[~,distanceErrorEntrancePupil,inputRayAngles] = evalCandidateLineOfSight(fixationPose,sceneGeometry,fixTargetWorld);

% Calculate and save the outputRay and the raypath
X = sceneGeometry.eye.landmarks.fovea.coords';
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;
[outputRay,rayPath] = rayTraceQuadrics(quadric.anglesToRay(X,rad2deg(inputRayAngles(1)),rad2deg(inputRayAngles(2))), opticalSystem);



end


%% LOCAL FUNCTIONS



%% calcDistanceFromEyeWorldTarget
% Performs the ray trace through the optical system of the
% eye and then calculates the Euclidean distance between the intersection
% point on the last surface a passed target in eyeWorld coordinates
function distanceError = calcDistanceFromEyeWorldTarget(opticalSystem,inputRay,eyeWorldTarget)
exitRay = rayTraceQuadrics(inputRay, opticalSystem);
distanceError = sqrt(sum((exitRay(:,1)-eyeWorldTarget).^2));
end


%% evalCandidateLineOfSight
% Given an eyePose, this routine finds the ray that departs the fovea and
% passes through the center of the entrance pupil, and then calculates how
% closely that ray passes to the fixationTarget in world coordinates.
function [distanceErrorFixationTarget, distanceErrorEntrancePupil, inputRayAngles] = evalCandidateLineOfSight(eyePose,sceneGeometry,fixTargetWorld)

% The "camera" of the sceneGeometry serves as the fixation target
sceneGeometry.cameraPosition.translation = fixTargetWorld;

% Obtain the center of the entrance pupil for this eye pose
[~, ~, ~, ~, eyePoints, pointLabels] = pupilProjection_fwd(eyePose, sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:));

% Define some options for the fmincon call in the loop
options = optimoptions(@fmincon,...
    'Display','iter');

% We will now find the initial angles of a ray leaving the fovea such that
% it intersects the center of the entrance pupil

% Define an error function that reflects the distance of the ray after it
% intersects the final surface from the center of the entrance pupil
X = sceneGeometry.eye.landmarks.fovea.coords';
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;
myError = @(p) calcDistanceFromEyeWorldTarget(opticalSystem,quadric.anglesToRay(X,p(1),p(2)),entrancePupilCenter');

% Supply an x0 guess which is the ray that connects the fovea with
% the center of the entrance pupil.
[~, angle_p1p2, angle_p1p3] = quadric.angleRays( [0 0 0; 1 0 0]', quadric.normalizeRay([X'; entrancePupilCenter-X']') );
angle_p1p3 = -angle_p1p3;
p0 = [angle_p1p2 angle_p1p3];

% Perform the search. This gets us the initial angles for a ray leaving the
% fovea and hitting the center of the entrance pupil.
[inputRayAngles, distanceErrorEntrancePupil] = fmincon(myError,p0,[],[],[],[],[-180,-180],[180,180],[],options);

% Now get the distance at which this ray intersects the fixation target in
% world coordinates. The local function that is called is copied from the
% virtualImageFunc routine.
distanceErrorFixationTarget = calcTargetIntersectError(entrancePupilCenter, ...
    inputRayAngles(1), inputRayAngles(2), eyePose, fixTargetWorld, ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.refraction.stopToCamera.opticalSystem);
end




%% calcTargetIntersectError
function distance = calcTargetIntersectError(eyePoint, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystem)
% Smallest distance of the exit ray from worldTarget
%
% Syntax:
%  distance = calcTargetIntersectError(eyePoint, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystem)
%
% Description:
%   This function returns the Euclidean distance between a target in the
%   world coordinate system and a ray that has exited from the optical
%   system of a rotated eye. This distance is calculated within eyeWorld
%   coordinates.
%
%   This function is used to find angles in the p1p2 and p1p3 planes that
%   minimize the distance between the intersection point of the ray in the
%   camera plane and the pinhole aperture of a camera. At a distance of
%   zero, the ray would enter the pinhole aperture of the camera and thus
%   produce a point on the resulting image.
%
% Inputs:
%   eyePoint
%   angle_p1p2, angle_p1p3 - Scalar in radians. The angle w.r.t. the 
%                           optical axis of the initial ray. 
%   eyePose               - As defined in the main function.
%   worldTarget
%   rotationCenters
%   opticalSystem
%
% Outputs:
%   distance              - Scalar in units of mm. The minimum Euclidean 
%                           distance between the worldTarget and a ray
%                           exiting the optical system of the rotated eye.
%                           Set to Inf if an error is returned by
%                           rayTraceQuadrics.
%

% Assemble the input ray. Note that the rayTraceQuadrics routine handles
% vectors as a 3x2 matrix, as opposed to a 2x3 matrix in this function.
% Tranpose operations ahead.
inputRay = quadric.anglesToRay(eyePoint',rad2deg(angle_p1p2),rad2deg(angle_p1p3));

% Ray trace
outputRayEyeWorld = rayTraceQuadrics(inputRay, opticalSystem);
outputRayEyeWorld = outputRayEyeWorld';

% If any must intersect surfaces were missed, the output ray will contain
% nans. In this case, return Inf for the distance
if any(isnan(outputRayEyeWorld))
    distance = Inf;
    return
end

% We assign the variable ET the coordinates of the target after conversion
% to eye coordinates. Then, the point is counter-rotated by the eye pose,
% so that the ET is in a position equivalent to if the eye had rotated.
cameraRot = -eyePose;
RotAzi = [cosd(cameraRot(1)) -sind(cameraRot(1)) 0; sind(cameraRot(1)) cosd(cameraRot(1)) 0; 0 0 1];
RotEle = [cosd(-cameraRot(2)) 0 sind(-cameraRot(2)); 0 1 0; -sind(-cameraRot(2)) 0 cosd(-cameraRot(2))];
RotTor = [1 0 0; 0 cosd(cameraRot(3)) -sind(cameraRot(3)); 0 sind(cameraRot(3)) cosd(cameraRot(3))];

% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. This is now the eyeTarget or ET
ET = worldTarget([3 1 2])';

% Torsion
ET=ET-rotationCenters.tor;
ET = (RotTor*ET')';
ET=ET+rotationCenters.tor;
% Elevation
ET=ET-rotationCenters.ele;
ET = (RotEle*ET')';
ET=ET+rotationCenters.ele;
% Azimuth
ET=ET-rotationCenters.azi;
ET = (RotAzi*ET')';
ET=ET+rotationCenters.azi;

% Calculate the distance between the closest approach of the outputRay to
% the target.
d = norm(cross(outputRayEyeWorld(2,:),ET - outputRayEyeWorld(1,:))) ...
    / norm(outputRayEyeWorld(2,:));
      
% Obtain the Euclidean distance in the 3 dimensions.
distance = sqrt(sum(d.^2));

end % calcTargetIntersectError

