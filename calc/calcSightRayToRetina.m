function [rayPath,retinaDistanceError] = calcSightRayToRetina(eye,rayDestination,rayOriginDistance,stopRadius,cameraMedium,surfaceTol)
% Ray that passes through entrance pupil center to reach the retinal coord
%
% Syntax:
%  [rayPath,retinaDistanceError] = calcSightRayToRetina(eye,rayDestination,rayOriginDistance,stopRadius,cameraMedium)
%
% Description
%   Given an eye and a coordinate on the retina, the routine identifies the
%   ray that passes through the center of the entrance window, and arrives
%   at the specified retinal location. If the retinal location is the
%   fovea, then this ray is the "line of sight" for the eye.
%
%   Note on terminology: the "entrance pupil" is the appearance of the
%   aperture stop through the optical elements that are in front of the
%   stop (e.g., the cornea of the eye) as viewed from a point on the
%   optical axis in object spaace. When the viewing point is off the
%   optical axis, then the image of the aperture stop is the "entrance
%   window".
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
%                           the ray from the incident node. Assumed to be
%                           1500 mm if not defined.
%   stopRadius            - Scalar. Radius of the aperture stop, in mm.
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
    [rayPath,retinaDistanceError] = calcSightRayToRetina(eye,rayDestination);
    % Confirm that the elements of the error vector are within tolerance
    assert(retinaDistanceError(1)<1e-3)
%}


arguments
    eye (1,1) {isstruct}
    rayDestination (3,1) {mustBeNumeric}
    rayOriginDistance (1,1)  {mustBeNumeric} = 1500
    stopRadius (1,1) {mustBeNumeric} = 1.53
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
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Check that the candidate beta and omega values are in range
    if abs(rayDestination(1))>90 || abs(rayDestination(2))>180
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Looks like a valid geodetic coordinate. Convert to Cartesian.
    rayDestination = quadric.ellipsoidalGeoToCart( rayDestination, S );
        
end

% Create a default sceneGeometry
sceneGeometry = createSceneGeometry('eye',eye,...
    'cameraMedium',cameraMedium,...
    'surfaceSetName',{'stopToMedium','mediumToCamera','mediumToRetina'});

% Define an objective, which is the distance of the point of intersection
% of the ray in the retina from the target coordinate.
myObj = @(p) objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius);

% p0 is set to direct from the emergent node to the retinal coordinate
referenceCoord = eye.landmarks.emergentNode.coords';
[p0(1),p0(2)] = quadric.rayToAngles(quadric.normalizeRay([rayDestination,referenceCoord-rayDestination]));

% Bounds
lb = [-90,-90];
ub = [ 90, 90];

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search. These angles identify a rayOrigin location which is at the
% desired angles with respect to the incident node.
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective with the final origin location
[retinaDistanceError,rayPath] = ...
    objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius);


end


%% Local function
function [fVal,rayPath] = objective(p,sceneGeometry,rayDestination,rayOriginDistance,stopRadius)

% The passed p vector is interpreted as the angles of a coordinate point
% w.r.t. the incident node. We find the coordinate at this point, at the
% rayOriginDistance
referenceCoord = sceneGeometry.eye.landmarks.incidentNode.coords';
R = quadric.anglesToRay(referenceCoord,p(1),p(2));
rayOrigin = R(:,1)+R(:,2).*rayOriginDistance;

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
