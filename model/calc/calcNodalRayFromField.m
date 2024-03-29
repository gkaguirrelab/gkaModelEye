function [rayPath,angleError] = calcNodalRayFromField(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,cameraMedium)
% Returns a nodal ray that arises from the specified field location
%
% Syntax:
%  [rayPath,angleError] = calcNodalRayFromField(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,cameraMedium)
%
% Description
%   Given an optical system (or eye structure) and a field location, the
%   routine returns a matrix that contains the path of a ray for which the
%   incident and emergent segments are parallel. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The fieldAngularPosition (and rayOriginDistance) is defined w.r.t. the
%   coordinate specified in angleReferenceCoord and distanceReferenceCoord.
%   If not defined, this is the origin of the longitudinal axis. Visual
%   angle is often defined w.r.t. the position of the approximate incident
%   nodal point, this coordinate should be provided if the
%   fieldAngularPosition is to be interpreted as visual angle. The distance
%   of a stimulus is relevant to the accommodative state of the eye. As
%   accommodation distance is measured relative to the principal point of
%   the eye, this value may be used for the distanceReferenceCoord.
%
% Inputs:
%   opticalSystem         - Either an eye structure (from which a
%                           "mediumToRetina" optical system in air will be
%                           derived), or an mx19 matrix, where m is set by
%                           the key value opticalSystemNumRows. Each row
%                           contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical] degrees with respect to the
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
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   angleError            - Scalar. The departure from parallel of the 
%                           incident and emergent rays (deg)
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Pick a field location
    F = [10,20];
    % Find the nodal ray
    [rayPath,angleError] = calcNodalRayFromField(eye,F);
    % Confirm that the angleError is within tolerance
    assert(angleError<1e-3)
%}


arguments
    opticalSystem {mustBeOpticalSystemCapable}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0,0]
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    cameraMedium = 'air'
end


% Check if we have a compiled version of findNodalRay
if exist('findNodalRayMex','file')==3
    findNodeHandle = @findNodalRayMex;
else
    findNodeHandle = @findNodalRay;
end

% Create the optical system
opticalSystem = parseOpticalSystemArgument(opticalSystem,'mediumToRetina',cameraMedium);

% Ensure that the opticalSystem has 100 rows to make the MEX file
% happy
osRowLength = size(opticalSystem,2);
opticalSystem = [opticalSystem; ...
    nan(100-size(opticalSystem,1),osRowLength)];

% Define the rayOrigin
fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);
rayOrigin = fieldRay(:,1);

% Find the nodal ray
[rayPath,~,angleError] = findNodeHandle(rayOrigin',opticalSystem);


end