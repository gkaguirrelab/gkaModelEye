function [rayPath,nodalPoints,errors] = calcSightRayToRetina(eye,rayDestination,rayOriginDistance,stopRadius,cameraMedium)
% Ray that passes through entrance pupil center to reach the retinal coord
%
% Syntax:
%  [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcSightRayToRetina(sceneGeometry,stopRadius,fixTargetDistance)
%
% Description
%   Given an eye and a retinal coordinate, the routine identifies the ray
%   that passes through the center of the entrance pupil, and arrives at
%   the specified retinal location.
%
%   If the retinal location is the fovea, then this ray is the "line of
%   sight" for the eye.
%
%   If not defined, the radius of the aperture stop is set to provide an
%   entrance pupil diameter of ~3.5 mm. The ray origin distance is assumed
%   to 500 mm unless set.
%
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
%   stopRadius            - Scalar. Radius of the aperture stop, in mm.
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
%   errors                - 1x2 matrix with the follow error values:
%                             - L2 norm of the distance between the retinal
%                               target coordinate and the intersection of
%                               the ray upon the retina
%                             - L2 norm of the mismatch between the desired
%                               and obtained rayOriginDistance.
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Obtain the coordinates of the fovea
    rayDestination = eye.landmarks.fovea.coords;
    % Find the sight ray to the fovea (i.e., the line of sight axis)
    [rayPath,nodalPoints,errors] = calcSightRayToRetina(eye,rayDestination);
    % Confirm that the first elements of the error vector are within
    % tolerance.
    assert(all(errors<1e-3))
%}


% Parse inputs
if nargin<2
    error('calcSightRayToRetina:invalidArguments','Too few input arguments')
end

if nargin==2
    rayOriginDistance = 1500;
    stopRadius = 1.53;
    cameraMedium = 'air';
end

if nargin==3
    stopRadius = 1.53;
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
end

% Make rayDestination a column vector
if all(size(rayDestination)==[1 3])
    rayDestination = rayDestination';
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
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Check that the candidate beta and omega values are in range
    if abs(rayDestination(1))>90 || abs(rayDestination(2))>180
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
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

% Create a default sceneGeometry (which has a pinhole camera)
sceneGeometry = createSceneGeometry('eye',eye,...
    'cameraMedium',cameraMedium,...
    'surfaceSetName',{'stopToMedium','mediumToCamera','mediumToRetina'});

% Define an objective, which is the distance of the point of intersection
% of the ray in the retina from the target coordinate.
myObj = @(p) objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius,findNodeHandle);

% p0 is set to direct from the retinal coordinate, through the location of
% an un-rotated corneal apex
[p0(1),p0(2)] = quadric.rayToAngles(quadric.normalizeRay([rayDestination,[0;0;0]-rayDestination]));

% Bounds
lb = [-90,-90];
ub = [ 90, 90];

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search. The angles that are returned from "p" are the angles of a ray
% w.r.t. the un-rotated corneal apex. If the search was successful, these
% angles identify a rayOrigin location which is at the desired angles with
% respect to the incident node.
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective with the final origin location
[retinalDistanceError,rayPath,nodalPoints,errors] = ...
    objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius,findNodeHandle);

% Consolidate the errors
errors = [retinalDistanceError, errors];

end


%% Local function
function [fVal,rayPath,nodalPoints,errors] = objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius,findNodeHandle)

% The passed p vector is interpreted as the angles of a coordinate point
% w.r.t. the origin of the coordinate system. We find the coordinate at
% this point, at the rayOriginDistance
rayOrigin = quadric.anglesToRay([0;0;0],p(1),p(2)).*rayOriginDistance;
rayOrigin = rayOrigin(:,2);

% Find the nodal points from this rayOrigin
[~,nodalPoints] = findNodeHandle(rayOrigin',sceneGeometry.refraction.mediumToRetina.opticalSystem);

% Adjust the rayOrigin so that it is at the appropriate rayOriginDistance
% w.r.t. the incident node.
rayOrigin = (rayOrigin-nodalPoints(:,1)).*(rayOriginDistance/norm(rayOrigin-nodalPoints(:,1)))+nodalPoints(:,1);

% Find the nodal points from this updated rayOrigin
[~,nodalPoints] = findNodeHandle(rayOrigin',sceneGeometry.refraction.mediumToRetina.opticalSystem);

% Obtain the L2 norm between the desired and realized rayOrigin
% distance, and store this in errors
errors = norm(norm(rayOrigin-nodalPoints(:,1))-rayOriginDistance);

% Obtain the center of the entrance pupil as seen from the rayOrigin
sceneGeometry.cameraPosition.translation = rayOrigin([2 3 1]);
[~, ~, ~, ~, ~, eyePoints, pointLabels] = projectModelEye([0 0 0 stopRadius], sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:))';

% Trace a ray from rayOrigin, intersecting the entrancePupilCenter, and
% follow it to the retina
R = quadric.normalizeRay([rayOrigin,entrancePupilCenter-rayOrigin]);
[outputRay, rayPath] = rayTraceQuadrics(R, sceneGeometry.refraction.mediumToRetina.opticalSystem);

% The fVal is the Euclidean distance between the retinal target coordinate
% and the intersection of the ray upon the retina.
fVal = norm(outputRay(:,1)-rayDestination);

end
