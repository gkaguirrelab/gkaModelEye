function [rayPath,nodalPoints,errors] = calcNodalRayToRetina(eye,rayDestination,rayOriginDistance,incidentNodeX0,cameraMedium)
% Returns the path of the nodal ray that intersects the retinal coordinate
%
% Syntax:
%  [rayPath,nodalPoints,errors] = calcNodalRayToRetina(eye,rayDestination,rayOriginDistance,incidentNodeX0,cameraMedium)
%
% Description
%   Given an eye structure and a retinal coordinate, the routine returns a
%   matrix that contains the path of a ray that arrives at this coord and
%   has an angle of incidence at cornea (w.r.t the optical axis) equal to
%   the angle with which it intersects the retina. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The nodal ray that arrives at the fovea is termed the "visual axis" of
%   the eye.
%
%   The routine returns the nodalPoints,which are found by extending the
%   initial and exit segments of the ray to their closest approach to the
%   optical axis. For a model eye with decentered and/or astigmatic
%   elements, there is not a single nodal point.
%
%   The routine can accept points on the ellipsoidal surface of the retina
%   specified in either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   rayDestination        - 1x3 or 3x1 vector that provides the coordinates
%                           of a location on the retina. The location can
%                           be given in Cartesian or geodetic coordinates.
%                           For the latter, the values are beta, omega, and
%                           elevation in units of degrees. Beta is defined
%                           over the range -90:90, and omega over the range
%                           -180:180. Elevation has an obligatory value of
%                           zero as this solution is only defined on the
%                           surface. The routine will attempt to determine
%                           which type of coordinate set has been provided.
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
%                             - L2 norm distance of ray intersection from 
%                               retinal target (in mm)
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
    % Obtain the coordinates of the fovea
    rayDestination = eye.landmarks.fovea.coords;
    % Find the nodal ray
    [rayPath,nodalPoints,errors] = calcNodalRayToRetina(eye,rayDestination);
    % Confirm that the first three elements of the error vector are within
    % tolerance. The final two elements are expected to be non-zero as a 
    % consequence of astigmatic and decentered elements in the model.
    assert(all(errors(1:3)<1e-2))
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

% Make  incidentNodeX0 and rayDestination column vectors
if all(size(rayDestination)==[1 3])
    rayDestination = rayDestination';
end
if all(size(incidentNodeX0)==[1 3])
    incidentNodeX0 = incidentNodeX0';
end

% Time to interpret the rayDestination variable. Hard-code a tolerance for
% the distance of the rayDestination coordinate from the retinal surface
surfaceTol = 1e-6;

% Get the retinal surface and quadric function
S = eye.retina.S;
Sfunc = quadric.vecToFunc(S);

% Interpret the value as Cartesian, if it is not a valid retinal location,
% the quadric function will evaluated with a non-zero value
if Sfunc(rayDestination(1),rayDestination(2),rayDestination(3)) > surfaceTol
    
    % If the last value of the coordinate is not zero, then this can't be a
    % geodetic coordinate either
    if rayDestination(3) > surfaceTol
        error('calcNodalRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Check that the candidate beta and omega values are in range
    if abs(rayDestination(1))>90 || abs(rayDestination(2))>180
        error('calcNodalRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Looks like a valid geodetic coordinate. Convert to Cartesian.
    rayDestination = quadric.ellipsoidalGeoToCart( rayDestination, S );
        
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

% Initialize an anonymous function for the objective.
myObj = @(p) objective(p,opticalSystem,rayDestination,rayOriginDistance,findNodeHandle);

% p0 is set to direct from the retinal coordinate, through a typical
% location for the incident node.
[p0(1),p0(2)] = quadric.rayToAngles(quadric.normalizeRay([rayDestination,incidentNodeX0-rayDestination]));

% Bounds
lb = [-90,-90];
ub = [ 90, 90];

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective function once more, using the final values
[retinalDistanceError,rayPath,nodalPoints,errors] = ...
    objective(p,opticalSystem,rayDestination,rayOriginDistance,findNodeHandle);

% Assemble the errors
errors = [retinalDistanceError,errors];


end


%% Local function

function [fVal,rayPath,nodalPoints,errors] = objective(p,opticalSystem,rayDestination,rayOriginDistance,findNodeHandle)
                                            
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

% Obtain the L2 norm between the desired and realized rayOrigin
% distance, and add this to the front of the errors
errors = [norm(norm(rayOrigin-nodalPoints(:,1))-rayOriginDistance), errors];

% Find the Euclidean distance between the retinal target coordinate and the
% intersection of the ray upon the retina.
fVal = norm(rayPath(:,end)-rayDestination);

end