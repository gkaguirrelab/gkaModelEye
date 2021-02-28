function [rayPath,angleError] = calcNodalRayFromField(opticalSystem,opticalFieldOrigin,rayOriginDistance,cameraMedium)
% Nodal ray that arises from the specified field location
%
% Syntax:
%  [rayPath,angleError] = calcNodalRayFromField(opticalSystem,opticalFieldOrigin,rayOriginDistance,cameraMedium)
%
% Description
%   Given an optical system (or eye structure) and field location, the
%   routine returns a matrix that contains the path of a ray for which the
%   incident and emergent segments are parallel. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   Note that opticalFieldOrigin is defined w.r.t. the origin of the
%   longitudinal axis. Visual angle is often defined w.r.t. the position of
%   the incident node (which itself is not a single point). Therefore, the
%   opticalFieldOrigin value specified here is not the same as the visual
%   angle for the same point in space.
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
%   opticalFieldOrigin    - 1x2 or 2x1 vector that provides the coordinates
%                           in degrees of the origin of the nodal ray with
%                           respect to the origin of the longitudinal axis
%                           (the un-rotated corneal apex).
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the longitudinal axis origin.
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


% Handle missing inputs
if nargin<2
    error('calcNodalRayFromField:invalidArguments','Too few input arguments')
end

if nargin==2
    rayOriginDistance = 1500;
    cameraMedium = 'air';
end

if nargin==3
    cameraMedium = 'air';
end

% Make opticalFieldOrigin a row vector
if all(size(opticalFieldOrigin)==[2 1])
    opticalFieldOrigin = opticalFieldOrigin';
end

% Check if we have a compiled version of findNodalRay
if exist('findNodalRayMex','file')==3
    findNodeHandle = @findNodalRayMex;
else
    findNodeHandle = @findNodalRay;
end

% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);
    end
end

% Define the rayOrigin
rayOrigin = quadric.anglesToRay([0;0;0],opticalFieldOrigin(1),opticalFieldOrigin(2)).*rayOriginDistance;
rayOrigin = rayOrigin(:,2);

% Find the nodal ray
[rayPath,angleError] = findNodeHandle(rayOrigin',opticalSystem);


end