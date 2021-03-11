function [rayPath,fieldAngularPosition,errors] = calcNodalRayToRetina(eye,rayDestination,rayOriginDistance,distanceReferenceCoord,incidentNodeX0,cameraMedium,surfaceTol)
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
%   The routine can accept points on the ellipsoidal surface of the retina
%   specified in either Cartesian or ellipsoidal geodetic coordinates.
%
%   The routine allows specification of the rayOriginDistance. This is of
%   little consequence for the path of the ray through the optics of the
%   eye, but is included here so that the calling function can define the
%   field origin of the rays with accuracy.
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
%                           1500 mm if not defined.
%   distanceReferenceCoord - 3x1 vector that provides the coordinate from
%                           which the rayOriginDistance is calculated. The
%                           The principal point is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   incidentNodeX0        - 3x1 vector that gives the location in eye
%                           space that is an initial guess for the
%                           location of the incident node of the optical
%                           system. If not supplied, a value that is
%                           typical for the human eye is used.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%   surfaceTol            - Scalar. When interpreting if the rayDestination
%                           is Cartesian or geodetic coordinates, this is
%                           the tolerance within which a candidate
%                           Cartesian coordinate must be on the quadric
%                           surface.
%
% Outputs:
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   fieldAngularPosition  - 1x2 vector with the angles of the origin of the
%                           nodal ray, w.r.t. the origin of the
%                           longitudinal axis [0;0;0]
%   errors                - 1x2 matrix with the follow error values:
%                             - L2 norm distance of ray intersection from 
%                               retinal target (in mm)
%                             - departure from parallel of the incident and
%                               emergent rays (deg)
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Obtain the coordinates of the fovea
    rayDestination = eye.landmarks.fovea.coords;
    % Find the visual axis
    [rayPath,fieldAngularPosition,errors] = calcNodalRayToRetina(eye,rayDestination);
    % Confirm that the errors are within tolerance.
    assert(all(errors<1e-2))
%}


arguments
    eye (1,1) {isstruct}
    rayDestination (3,1) {mustBeNumeric}
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0; 0; 0]
    incidentNodeX0 (3,1) {mustBeNumeric} = [-7; 0; 0]
    cameraMedium = 'air'
    surfaceTol (1,1) {mustBeNumeric} = 1e-6
end

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

% Create the optical system
opticalSystem = parseOpticalSystemArgument(eye,'mediumToRetina',cameraMedium);

% Initialize an anonymous function for the objective.
myObj = @(p) objective(p,opticalSystem,rayDestination,rayOriginDistance,distanceReferenceCoord,findNodeHandle);

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
fieldAngularPosition = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective function once more, using the final values
[retinalDistanceError,rayPath,angleError] = ...
    objective(fieldAngularPosition,opticalSystem,rayDestination,rayOriginDistance,distanceReferenceCoord,findNodeHandle);

% Assemble the errors
errors = [retinalDistanceError,angleError];


end


%% Local function

function [fVal,rayPath,angleError] = objective(fieldAngularPosition,opticalSystem,rayDestination,rayOriginDistance,distanceReferenceCoord,findNodeHandle)
                                            
% The passed p vector is interpreted as the angles of a coordinate point
% w.r.t. the origin of the coordinate system.
fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,[0;0;0],distanceReferenceCoord);
rayOrigin = fieldRay(:,1);

% Find the nodal ray from this point
[rayPath,~,angleError] = findNodeHandle(rayOrigin',opticalSystem);

% Find the Euclidean distance between the retinal target coordinate and the
% intersection of the ray upon the retina.
fVal = norm(rayPath(:,end)-rayDestination);

end