function chiefRayPath = calcChiefRay(sceneGeometry,eyePose)
% Returns path for the chief ray from the camera to the model eye
%
% Syntax:
%  chiefRayPath = calcChiefRay(sceneGeometry,eyePose)
%
% Description
%   The chief ray for a point in object space is the ray that travels from
%   that point through the center of the aperture stop of the optical
%   system.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The stop
%                           radius value is unused by this routine.
%
% Outputs:
%   chiefRayPath          - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for
%                           chiefRayPath(1,:) is equal to initial position.
%                           If a surface is missed, then the coordinates
%                           for that surface will be nan.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    eyePose = [0 0 0 2];
    chiefRayPath = calcChiefRay(sceneGeometry,eyePose);
%}

% Define the optical system
opticalSystem = sceneGeometry.refraction.stopToCamera.opticalSystem;

% Obtain the chief ray. This is a ray that originates at the stop center
% and arrives at the pinhole aperture of the camera. The output of this
% function is the chief ray as it intersects the stop.
[~,chiefRayFromStop] = virtualImageFunc( sceneGeometry.eye.stop.center, ...
    eyePose, sceneGeometry.cameraPosition.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    opticalSystem );

% Obtain the final segment of the ray. This is the ray as it departs the
% last surface of the system and is directed towards the aperture of the
% camera.
[~,chiefRayPath] = rayTraceQuadrics(chiefRayFromStop, opticalSystem);

end

