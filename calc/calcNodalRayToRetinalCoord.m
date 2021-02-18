function [rayPath,nodalPoints,errors] = calcNodalRayToRetinalCoord(eye,G,X,originRayDepth,cameraMedium)
% Returns the path of the nodal ray that intersects the retinal coordinate
%
% Syntax:
%  [rayPath,nodes,errors] = calcNodalRay(eye,G,X,originRayDepth,cameraMedium)
%
% Description
%   Given an eye structure and a coordinate, the routine returns a matrix
%   that contains the path of a ray that arrives at this coord and has an
%   angle of incidence at cornea (w.r.t the optical axis) equal to the
%   angle with which it intersects the retina. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The routine returns the nodalPoints,which are found by extending the
%   initial and exit segments of the ray to the optical axis. The
%   nodal ray that arrives at the fovea is termed the "visual axis" of the
%   eye.
%
%   For a model eye with decentered and/or astigmatic elements, there may
%   not exist a true nodal ray, but only approximations to this concept.
%   There may not be a ray that both has incident and emergent rays, and
%   intersects the retinal target. Further, the incident and emergent rays,
%   when extended, will pass by but not intersect the optical axis. These
%   various error are supplied in the errors output variable. Further, the
%   particular nodal points will vary depending upon the particular retinal
%   target.
%
%   The routine can accept points on the ellipsoidal surface of the retina
%   specified in either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   G                     - 1x3 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X                     - 1x3 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   originRayDepth        - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex.
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
%                             - distance of ray intersection from retinal
%                               target (in mm)
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
    % Pick a retinal point in geodetic coordinates, a bit away from the
    % vertex
    G = [-65,-65,0];
    % Find the nodal ray
    [rayPath,nodalPoints,errors] = calcNodalRayToRetinalCoord(eye,G);
    % Show the optical system, nodal ray, and nodal points
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'rayPath',rayPath,'surfaceAlpha',0.05);
    hold on
    xlim([-25 10])
    plot3(nodalPoints(1,:),nodalPoints(2,:),nodalPoints(3,:),'*b')
%}


% Parse inputs
if nargin<2
    error('Invalid number of input arguments');
end

if nargin==2
    X = [];
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==3
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
end

% Check if the X value is empty. If so, derive the X
% Cartesian coordinates from the ellipsoidal geodetic coordinate.
if isempty(X)
    S = eye.retina.S;
    X = quadric.ellipsoidalGeoToCart( G', S );
end

% Convert X to a column vector
if all(size(X)==[1 3])
    X = X';
end

% Check if we have a compiled version of findNodalRay
if exist('findNodalRayMex')==3
    findeNodeHandle = @findNodalRayMex;
else
    findeNodeHandle = @findNodalRay;
end

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);


% Initialize an anonymous function for the objective. This combines
% both the objective and the non-linear constraint.
myObj = @(x) objective(x,opticalSystem,originRayDepth,X,findeNodeHandle);


% x0 is placed so that it is in line with the retinal
% target location and the corneal apex, with the two values being the
% horizontal and vertical position of the ray origin (in mm)
T = (originRayDepth./X(1)) .* X;
T = [sqrt(originRayDepth^2-T(2)^2-T(3)^2);T(2);T(3)];
x0 = [T(2),T(3)];

% Search
x = fminsearch(myObj,x0);

% Evaluate the objective function once more, using the final values

[retinalDistanceError,rayPath,nodalPoints,errors] = ...
    objective(x,opticalSystem,originRayDepth,X,findeNodeHandle);

% Assemble the errors
errors = [retinalDistanceError,errors];


end


%% Local function

function [fVal,rayPath,nodalPoints,errors] = objective(x,opticalSystem,originRayDepth,X,findeNodeHandle)

% Define a point in space
T = [sqrt(originRayDepth^2-x(1)^2-x(2)^2);x(1);x(2)];

% Find the nodal ray from this point
[rayPath,nodalPoints,errors] = findeNodeHandle(T',opticalSystem);

% Find the Euclidean distance between the retinal target coordinate and the
% intersection of the ray upon the retina.
fVal = norm(rayPath(:,end)-X);

end