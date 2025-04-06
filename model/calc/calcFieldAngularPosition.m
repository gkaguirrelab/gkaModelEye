function [fieldAngularPosition,targetDistance] = calcFieldAngularPosition(eye,targetWorldCoordinate)
% Visual field coordinates of world coordinate fixation point
%
% Syntax:
%  [fieldAngularPosition,targetDistance] = calcFieldAngularPosition(eye,targetWorldCoordinate)
%
% Description
%   Given a world coordinate target, provides the location in the visual
%   field in horizontal and vertical degrees w.r.t. the longitudinal axis
%   of the eye in primary position, and the distance of that target in mm
%   from the incident node of the eye.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   targetWorldCoordinate - 3x1 vector that provides the world coordinates
%                           of a target.
%
% Outputs:
%   fieldAngularPosition  - 1x2 vector that provides the coordinates in
%                           degrees of visual angle of the target
%                           relative to the longitudinal axis of the eye
%                           when it is aligned with the camera.
%   targetDistance        - Scalar. The distance (in mm) of the target from
%                           the incident node.
%
% Examples:
%{
    eye = modelEyeParameters();
    targetWorldCoordinate = [50; 0; 1000];
    [fieldAngularPosition,targetDistance] = calcFieldAngularPosition(eye,targetWorldCoordinate);
%}


arguments
    eye (1,1) {isstruct}
    targetWorldCoordinate (3,1) {mustBeNumeric} = [0; 0; 1500]
end

referenceCoord = eye.landmarks.incidentNode.coords;
pp = [referenceCoord', convertWorldToEyeCoord(targetWorldCoordinate)'];
R = quadric.coordsToRay(pp);
[fieldAngularPosition(1),fieldAngularPosition(2)] = quadric.rayToAngles( R );
targetDistance = norm(pp);

end

