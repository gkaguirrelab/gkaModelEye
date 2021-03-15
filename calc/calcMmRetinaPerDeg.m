function mmPerDeg = calcMmRetinaPerDeg(eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,deltaDegEuclidean,cameraMedium)
% Returns mm/degree of field angle at a field location
%
% Syntax:
%  mmPerDeg = calcMmRetinaPerDeg(eye,fieldOrigin,deltaDegEuclidean,cameraMedium)
%
% Description
%   Given an eye structure and the location of a point in the field,
%   returns the mm of retina per degree of field angle at the corresponding
%   retinal location.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in angleReferenceCoord.
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
%   deltaDegEuclidean     - Scalar. The measurement is made for a small
%                           displacement of visual angle specified by this
%                           variable.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   mmPerDeg              - Scalar. The mm of retina per degree of visual
%                           angle at the specified visual field location.
%
% Examples:
%{
    % mm of retina / deg at the fovea
    eye = modelEyeParameters();
    fieldAngularPosition = eye.landmarks.fovea.degField;
    rayOriginDistance = 1500;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    mmPerDeg = calcMmRetinaPerDeg(eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord);
    fprintf('%2.3f retinal mm per deg visual field at the fovea in the emmetropic eye.\n',mmPerDeg);
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0 0]
    rayOriginDistance {isscalar,mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) double = [0;0;0]
    deltaDegEuclidean {isscalar,mustBeNumeric} = 1
    cameraMedium = 'air'
end


% Set up the jitter in angles around the specified field location
deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];

% Trace the nodal rays from this field position
rayPath0 = calcNodalRayFromField(eye,fieldAngularPosition-deltaAngles./2,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);
rayPath1 = calcNodalRayFromField(eye,fieldAngularPosition+deltaAngles./2,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);

% Find the geodetic distance between the points on the retina
geodesic = quadric.geodesic(eye.retina.S,[],[],rayPath0(:,end),rayPath1(:,end));

% Calculate the mm per deg
mmPerDeg = geodesic / norm(deltaAngles);


end