function [outputRay,rayPath,angleError] = calcNodalRay(eye,coord,coordType,traceDirection,cameraMedium)
% Returns the path of the nodal ray that intersects the coordinate
%
% Syntax:
%  [outputRay,rayPath, angleError] = calcNodalRay(eye,G,X,cameraMedium)
%
% Description
%   Given an eye structure and a coordinate, the routine returns a matrix
%   that contains the path of a ray that satisfies the property that the
%   angle (w.r.t the optical axis) of the ray as it departs the retina is
%   equal to the angle with which it departs the cornea. This is a "nodal
%   ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The passed coordinate is treated either as the origin or destination of
%   the ray, dependning upon the "traceDirection" setting. Note that the
%   nodal ray that arrives at the fovea is termed the "visual axis" of the
%   eye.
%
%   The routine can accept points on the ellipsoidal surface of the retina
%   specified in either Cartesian or ellipsoidal geodetic coordinates, or
%   from locations in the visual field specified in degrees or mm. The
%   interpretation of the coordinate is provided by the "coordType"
%   variable.
%
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   coord                 - 2x1 or 3x1 vector that provides the coordinates
%                           of the origin point of the ray to be traced.
%                           This coordinate can be either within the eye or
%                           in the visual field, and can be of the form:
%                             - retinal geodetic (beta, omega, and
%                               optionally elevation in units of degrees.
%                               Beta is defined over the range -90:90, and
%                               omega over the range -180:180. Elevation
%                               has an obligatory value of zero as this
%                               solution is only defined on the surface.
%                             - visual field (horizontal, vertical) in
%                               units of degrees.
%                             - cartesian (x,y,z) in units of mm.
%	coordType             - Char vector or string. Defines the nature of
%                           the coord variable. Options are:
%                               {'geodetic','field','cartesian'}
%   traceDirection        - Char vector or string. Options are:
%                               {'cameraToEye','eyeToCamera'}
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   outputRay             - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   angleError            - Scalar. The angle between the initial and
%                           output rays for the nominal nodal ray. Ideally,
%                           this value should be zero.
%
% Examples:
%{
    eye = modelEyeParameters();    
    [outputRay,rayPath,angleError] = calcNodalRay(eye,[],eye.landmarks.vertex.coords);
    assert(angleError < 1e-3);
%}


% Handle incomplete arguments
if nargin<4
    error('Invalid number of input arguments');
end

if nargin==4
    cameraMedium = 'air';
end

% Convert geodetic and 


if nargin==2
    % If only two input values were passed, derive the X0 Cartesian
    % coordinates from the ellipsoidal geodetic coordinates.
    S = eye.retina.S;
    X = quadric.ellipsoidalGeoToCart( G, S );
    cameraMedium = 'air';
end

if nargin==3
    % Check if the X0 value is empty. If so, derive the X0
    % Cartesian coordinates from the ellipsoidal geodetic coordinate.
    if isempty(X)
        S = eye.retina.S;
        X = quadric.ellipsoidalGeoToCart( G, S );
    end
    cameraMedium = 'air';
end

if nargin==4
    % Check if the X0 value is empty. If so, derive the X0
    % Cartesian coordinates from the ellipsoidal geodetic coordinate.
    if isempty(X)
        S = eye.retina.S;
        X = quadric.ellipsoidalGeoToCart( G, S );
    end
end

% Make sure that X is a column vector
if all(size(X)==[1 3])
    X = X';
end

% Assemble the optical system
opticalSystem = assembleOpticalSystem( eye, ...
    'surfaceSetName','retinaToMedium', 'cameraMedium', cameraMedium);

% Define an error function that reflects the difference in angles from the
% initial ray and the output ray from the optical system
myError = @(p) calcOffsetFromParallel(opticalSystem,quadric.anglesToRay(X,p(1),p(2)));

% Supply an x0 guess. Check a few targets to try and find one that produces a valid
% trace through the optical system.
targets = [...
    eye.lens.back', ...                      % lens posterior
    ((eye.lens.back+eye.lens.center)./2)',...  % mid-point of posterior lens
    eye.lens.center'];                       % lens center

stillSearching = true;
targetIdx = 1;
while stillSearching
    R = quadric.normalizeRay([X'; (targets(:,targetIdx)-X)']');
    [outputRay,rayPath] = rayTraceQuadrics(R, opticalSystem);    
    if any(isnan(outputRay))
        targetIdx = targetIdx + 1;
        if targetIdx > length(targets)
            angleError = nan;
            return
        end
    else
        stillSearching = false;
    end
end
[angle_p1p2, angle_p1p3] = quadric.rayToAngles( R );
x0 = [angle_p1p2, angle_p1p3];

% Add a small offset to x0, as the search can be stuck in a local minimum
% if it starts pointed straight down the optical system
x0 = x0+[0.01 0.01];

% define some search options
options = optimset('Display','off');

% Perform the search
[inputRayAngles, angleError] = fminsearch(myError,x0, options);

% Calculate and save the outputRay and the raypath
[outputRay,rayPath] = rayTraceQuadrics(quadric.anglesToRay(X,inputRayAngles(1),inputRayAngles(2)), opticalSystem);

end


%% Local functions


% Local function. Performs the ray trace through the optical system of the
% eye and then calculates the angle between the initial ray and the output
% ray.
function angleError = calcOffsetFromParallel(opticalSystem,inputRay)
exitRay = rayTraceQuadrics(inputRay, opticalSystem);
angleError = quadric.angleRays( inputRay, exitRay );
end
