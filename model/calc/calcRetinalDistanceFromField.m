function [distance, geodesicPathCoords] = calcRetinalDistanceFromField(opticalSystem,fieldAngularPositionStart,fieldAngularPositionEnd,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium)
% Geodesic distance between retinal points defined by two field positions
%
% Syntax:
%  [distance, geodesicPathCoords] = calcRetinalDistanceFromField(opticalSystem,fieldAngularPositionStart,fieldAngularPositionEnd,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium)
%
% Description
%   Given an optical system (or eye structure) and two field locations, the
%   routine returns the geodesic distance between the two points on the
%   retina that subtend the field locations.
%
% Inputs:
%   opticalSystem         - A sceneGeometry, eye structure, optical system
%                           structure, or optical system matrix. If an eye
%                           or scene structure is supplied, then the
%                           "mediumToRetina" optical system in air will be
%                           derived).
%   fieldAngularPositionStart, fieldAngularPositionEnd - 2x1 vectors that 
%                           provide the coordinates of points in the visual
%                           field in [horizontal, vertical] degrees with
%                           respect to longitudinal axis of the optical
%                           system.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the distanceReferenceCoord.
%   angleReferenceCoord   - 3x1 vector that provides the coordinate from
%                           which the ray origin angles and distance are
%                           to be calculated. By default, this is [0;0;0],
%                           which is the origin coordinate on the
%                           longitudinal axis.
%   distanceReferenceCoord - 3x1 vector that provides the coordinate from
%                           which the rayOriginDistance is calculated. The
%                           The principal point is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   distance              - Scalar. Approximation to the geodesic distance   
%                           (in mm) between the two retinal points.
%   geodesicPathCoords    - 3x50 matrix of locations along the geodesic.
%
% Examples:
%{
    % Distance from the fovea to the optic disc on the retina
    eye = modelEyeParameters();
    F1 = eye.landmarks.fovea.degField;
    F2 = eye.landmarks.opticDisc.degField;
    d = calcRetinalDistanceFromField(eye,F1,F2);
%}


arguments
    opticalSystem {mustBeOpticalSystemCapable}
    fieldAngularPositionStart (1,2) {mustBeNumeric}
    fieldAngularPositionEnd (1,2) {mustBeNumeric} = [0,0]
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    cameraMedium = 'air'
end

% Create the optical system
opticalSystem = parseOpticalSystemArgument(opticalSystem,'mediumToRetina',cameraMedium);

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Find the retinal points
rayPathA = calcNodalRayFromField(opticalSystem,fieldAngularPositionStart,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);
rayPathB = calcNodalRayFromField(opticalSystem,fieldAngularPositionEnd,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);

% Obtain the geodesic distance on the last surface in the optical system
[distance, geodesicPathCoords] = quadric.geodesic(opticalSystem(end,1:10),[rayPathA(:,end),rayPathB(:,end)]);

end