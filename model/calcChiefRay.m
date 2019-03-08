function [outputRay,rayPath] = calcChiefRay(sceneGeometry,eyePose)
% Returns path for the chief ray from the camera to the model eye
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
    sceneGeometry = createSceneGeometry();
    eyePose = [0 0 0 2];
    [outputRay,rayPath]=calcChiefRay(sceneGeometry,eyePose);
%}

% Define the optical system
opticalSystem = sceneGeometry.refraction.stopToCamera.opticalSystem;

% Assemble the args for the virtualImageFunc
args = {sceneGeometry.cameraPosition.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    opticalSystem};

% Obtain the chief ray. This is a ray that originates at the stop
% center and arrives at the pinhole aperture of the camera. The output of
% this function is the chief ray as it intersects the stop.
[~,chiefRayFromStop] = virtualImageFunc( sceneGeometry.eye.stop.center, [0 0 0 2], args{:} );

% Obtain the final segment of the ray. This is the ray as it departs the
% last surface of the system and is directed towards the aperture of the
% camera.
chiefRayFromLastSurface = rayTraceQuadrics(chiefRayFromStop, opticalSystem);

% Find the closest point to the optical axis of the
% extension of the chiefRayFromLastSurface. This point provides axial position
% along the optical axis at which 

[centerPoint, distance, R1p, R2p]=distanceRays(R1,R2)


end

