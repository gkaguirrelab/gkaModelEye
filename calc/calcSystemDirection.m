function systemDirection = calcSystemDirection(opticalSystem, rayStartDepth)
% Returns the valid direction for ray tracing for an optical system
%
% Syntax:
%  systemDirection = calcSystemDirection(opticalSystem)
%
% Description
%   The implementation of optical systems and ray tracing in this code
%   results in only one direction of ray tracing being available for a
%   given opticalSystem variable. We try both directions here and report
%   back the valid solution.
%
% Inputs:
%   opticalSystem         - An mx19 matrix, where m is set by the key value
%                           opticalSystemNumRows. Each row contains the
%                           values:
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
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex. Assumed to be
%                           1e5 mm if not defined.
%
% Outputs:
%  'systemDirection'      - Char vector with valid values 'eyeToCamera' or
%                           'cameraToEye'. Defines the direction of ray
%                           tracing for this optical system. For the
%                           special case of an optical system that is empty
%                           except for the specification of the refractive
%                           index of the medium, the routine returns
%                           'cameraToEye; eyeToCamera'.
%
% Examples:
%{
    % Test if a system created in a particular direction is correctly classed
    sceneGeometry = createSceneGeometry();
    systemDirectionIn = 'eyeToCamera';
    opticalSystem = sceneGeometry.refraction.stopToMedium.opticalSystem;
    systemDirectionOut = calcSystemDirection(opticalSystem);
    assert(strcmp(systemDirectionIn,systemDirectionOut));
%}

% Handle nargin
if nargin==1
    rayStartDepth = 1e5;
end

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% First check to make sure that the opticalSystem variable is well formed
if size(opticalSystem,2)~=19
    systemDirection = 'Not valid matrix dimensions';
    return
end
if ~all(isnan(opticalSystem(1,1:18)))
    systemDirection = 'Missing initial state row';
    return
end
if any(any(isnan(opticalSystem(2:end,:))))
    systemDirection = 'Invalid nan in matrix';
    return
end

% Handle the special case of an optical system that is empty except for
% specification of the refractive index of the medium
if size(opticalSystem,1)==1
    systemDirection = 'cameraToEye; eyeToCamera';
    return
end

% Trace an axial ray from the right (cameraToEye)
R1 = quadric.normalizeRay(quadric.anglesToRay([rayStartDepth;0;0],180,0));
M1 = rayTraceQuadrics(R1, opticalSystem);

% Trace an axial ray from the left (eyeToCamera)
R2 = quadric.normalizeRay(quadric.anglesToRay([-rayStartDepth;0;0],0,0));
M2 = rayTraceQuadrics(R2, opticalSystem);

% Find the non-nan value
if ~isnan(M1)
    systemDirection = 'cameraToEye';
elseif ~isnan(M2)
    systemDirection = 'eyeToCamera';
else
    systemDirection = 'No valid ray trace on optical axis';
end


end

