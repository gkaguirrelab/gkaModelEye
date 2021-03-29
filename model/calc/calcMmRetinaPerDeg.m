function mmPerDeg = calcMmRetinaPerDeg(eye,fieldAngularPosition,deltaAngles,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium)
% Returns mm/degree of field angle at a field location
%
% Syntax:
%  mmPerDeg = calcMmRetinaPerDeg(eye,fieldAngularPosition,deltaAngles,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium)
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
%   deltaAngles           - 1x2 vector. The measurement is made for a small
%                           displacement of horizontal and vertical visual
%                           angle specified by this variable. Setting this
%                           value allows one to measure the unit conversion
%                           in different polar directions around the
%                           specified point.
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
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   mmPerDeg              - Scalar. The mm of retina per degree of visual
%                           angle at the specified visual field location.
%
% Examples:
%{
    % Simple example for the Drasdo model eye
    eye = modelEyeParameters('species','Drasdo');
    mmPerDeg = calcMmRetinaPerDeg(eye);
%}
%{
    % mm of retina / deg at the fovea
    eye = modelEyeParameters();
    fieldAngularPosition = eye.landmarks.fovea.degField;
    deltaAngles = [1 0]; % Measure horizontally separated points
    rayOriginDistance = 1500;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    mmPerDeg = calcMmRetinaPerDeg(eye,fieldAngularPosition,deltaAngles,rayOriginDistance,angleReferenceCoord);
    fprintf('%2.3f retinal mm per deg visual field at the fovea in the emmetropic eye.\n',mmPerDeg);
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0 0]
    deltaAngles (1,2) {mustBeNumeric} = [1/sqrt(2) 1/sqrt(2)]
    rayOriginDistance {isscalar,mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) double = [0;0;0]
    cameraMedium = 'air'
end


% Trace the nodal rays from this field position
rayPath0 = calcNodalRayFromField(eye,fieldAngularPosition-deltaAngles/2,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);
rayPath1 = calcNodalRayFromField(eye,fieldAngularPosition+deltaAngles/2,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium);

% Assemble the points
P = [rayPath0(:,end),rayPath1(:,end)];

% Find the geodetic distance between the points on the retina
geodesic = quadric.geodesic(eye.retina.S,P);

% Calculate the mm per deg
mmPerDeg = geodesic / norm(deltaAngles);


end