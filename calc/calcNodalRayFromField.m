function [rayPath,nodalPoints,errors] = calcNodalRayFromField(eye,fieldOrigin,rayOriginDistance,incidentNodeX0,cameraMedium)
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
%                           the ray from the corneal apex. Assumed to be
%                           500 mm if not defined.
%   incidentNodeX0        - An optional 1x3 vector that gives the location
%                           in eye  space that is an initial guess for the
%                           location of the incident node of the optical
%                           system. If not supplied, a value that is
%                           typical for the human eye is used.
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
    % Show the optical system, nodal ray, and nodal points
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'rayPath',rayPath,'surfaceAlpha',0.05);
    hold on
    xlim([-25 10])
    plot3(nodalPoints(1,:),nodalPoints(2,:),nodalPoints(3,:),'*b')
%}


% Handle missing inputs
if nargin<2
    error('calcNodalRayFromField:invalidArguments','Too few input arguments')
end

if nargin==2
    rayOriginDistance = 500;
    incidentNodeX0 = [-7 0 0];
    cameraMedium = 'air';
end

if nargin==3
    incidentNodeX0 = [-7 0 0];
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
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

% An anonymous function that returns the coordinates of a point that is at
% the specified rayOriginDistance, and at the passed angles (w.r.t. to the
% un-rotated corneal apex)
my2ndColumn = @(R) R(:,2);
myCoord = @(p) my2ndColumn(quadric.anglesToRay([0;0;0],p(1),p(2)).*rayOriginDistance);

% Initialize an anonymous function for the objective.
myObj = @(p) objective(myCoord(p),opticalSystem,fieldOrigin,findNodeHandle);

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
    objective(myCoord(p),opticalSystem,fieldOrigin,findNodeHandle);

% Assemble the errors
errors = [angleMatchError,errors];


end


%% Local function

function [fVal,rayPath,nodalPoints,errors] = objective(rayOrigin,opticalSystem,fieldOrigin,findeNodeHandle)

% Find the nodal ray from this point
[rayPath,nodalPoints,errors] = findeNodeHandle(rayOrigin',opticalSystem);

% Find the angle of the rayOrigin with respect to the incident node
[vf(1), vf(2)] = quadric.rayToAngles(quadric.normalizeRay([nodalPoints(:,1),rayOrigin-nodalPoints(:,1)]));

% Obtain the L2 norm of the mis-match between the desired and obtained
% visual field location
fVal = norm(fieldOrigin-vf);

end