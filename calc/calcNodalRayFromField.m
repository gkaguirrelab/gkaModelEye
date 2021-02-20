function [rayPath,nodalPoints,errors] = calcNodalRayFromField(eye,fieldOrigin,rayOriginDistance,cameraMedium)
% Nodal ray that arises from the specified visual field location
%
% Syntax:
%  [rayPath,nodalPoints,errors] = calcNodalRayFromField(eye,fieldOrigin,rayOriginDistance,incidentNodeX0,cameraMedium)
%
% Description
%   Given an eye structure and visual field location, the routine returns a
%   matrix that contains the path of a ray that arises from this location
%   and has an angle of incidence at cornea (w.r.t the optical axis) equal
%   to the angle with which it intersects the retina. This is a "nodal
%   ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   Visual angle is typically defined as the angle of an location with
%   respect to the incident node of an optical system. For an astigmatic,
%   decentered optical system, there is not a single nodal point.
%   Consequently, this routine conducts an interative search to find a ray
%   that has the specified visual angle with respect to its own incident
%   node.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   fieldOrigin           - 1x2 or 2x1 vector that provides the coordinates
%                           in degrees of visual angle of the origin of the
%                           nodal ray.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the incident node. Assumed to be
%                           500 mm if not defined.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   nodalPoints           - 3x2 matrix that provides the approximation to
%                           incident and emergent nodal points in the eye
%                           coordinate space. This is the point on each ray
%                           that is closest to the optical axis.
%   errors                - 1x4 matrix with the follow error values:
%                             - L2 norm of the mismatch between the desired
%                               and obtained visual angles.
%                             - L2 norm of the mismatch between the desired
%                               and obtained rayOriginDistance.
%                             - departure from parallel of the incident and
%                               emergent rays (deg)
%                             - distance of the incident nodal point from
%                               the incident ray
%                             - distance of the emergent nodal point from
%                               the emergent ray
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Pick a visual field location
    F = [10,20];
    % Find the nodal ray
    [rayPath,nodalPoints,errors] = calcNodalRayFromField(eye,F);
    % Confirm that the first three elements of the error vector are within
    % tolerance. The final two elements are expected to be non-zero as a 
    % consequence of astigmatic and decentered elements in the model.
    assert(all(errors(1:3)<1e-3))
%}


% Handle missing inputs
if nargin<2
    error('calcNodalRayFromField:invalidArguments','Too few input arguments')
end

if nargin==2
    rayOriginDistance = 500;
    cameraMedium = 'air';
end

if nargin==3
    cameraMedium = 'air';
end

% If the length of fieldOrigin is 3, and the last element is zero, drop
% this as it is a torsion place holder.
if length(fieldOrigin)==3
    if fieldOrigin(end)==0
        fieldOrigin = fieldOrigin(1:2);
    else
        error('calcNodalRayFromField:invalidArguments','Field origin should be two elements')
    end
end

% Make fieldOrigin a row vector
if all(size(fieldOrigin)==[2 1])
    fieldOrigin = fieldOrigin';
end

% Check if we have a compiled version of findNodalRay
if exist('findNodalRayMex','file')==3
    findNodeHandle = @findNodalRayMex;
else
    findNodeHandle = @findNodalRay;
end

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Initialize an anonymous function for the objective. The fVal is the error
% in matching the desired fieldOrigin location.
myObj = @(p) objective(p,opticalSystem,fieldOrigin,rayOriginDistance,findNodeHandle);

% Bounds
lb = [-90,-90];
ub = [ 90, 90];

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search. The desired fieldOrigin itself is used as the p0 value.
% The angles that are returned from "p" are the angles of a ray w.r.t. the
% un-rotated corneal apex. If the search was successful, these angles
% identify a rayOrigin location which is at the desired angles with respect
% to the incident node.
p = fmincon(myObj,fieldOrigin,[],[],[],[],lb,ub,[],options);

% Evaluate the objective function once more, using the final values
[angleMatchError,rayPath,nodalPoints,errors] = ...
    objective(p,opticalSystem,fieldOrigin,rayOriginDistance,findNodeHandle);

% Assemble the errors
errors = [angleMatchError,errors];


end


%% Local function

function [fVal,rayPath,nodalPoints,errors] = objective(p,opticalSystem,fieldOrigin,rayOriginDistance,findNodeHandle)

% The passed p vector is interpreted as the angles of a coordinate point
% w.r.t. the origin of the coordinate system. We find the coordinate at
% this point, at the rayOriginDistance
rayOrigin = quadric.anglesToRay([0;0;0],p(1),p(2)).*rayOriginDistance;
rayOrigin = rayOrigin(:,2);

% Find the nodal points from this point
[~,nodalPoints] = findNodeHandle(rayOrigin',opticalSystem);

% Adjust the rayOrigin so that it is at the appropriate rayOriginDistance
% w.r.t. the incident node.
rayOrigin = (rayOrigin-nodalPoints(:,1)).*(rayOriginDistance/norm(rayOrigin-nodalPoints(:,1)))+nodalPoints(:,1);

% Repeat the ray trace from this updated point
[rayPath,nodalPoints,errors] = findNodeHandle(rayOrigin',opticalSystem);

% Find the angle of the rayOrigin with respect to the incident node
[vf(1), vf(2)] = quadric.rayToAngles(quadric.normalizeRay([nodalPoints(:,1),rayOrigin-nodalPoints(:,1)]));

% Obtain the L2 norm between the desired and realized rayOrigin
% distance, and add this to the front of the errors
errors = [norm(norm(rayOrigin-nodalPoints(:,1))-rayOriginDistance), errors];

% Obtain the L2 norm of the mis-match between the desired and obtained
% visual field location
fVal = norm(fieldOrigin-vf);

end