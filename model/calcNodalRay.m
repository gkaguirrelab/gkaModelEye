function [outputRay,rayPath, angleError] = calcNodalRay(eye,G,X,cameraMedium)
% Returns the path of the nodal ray from a retinal point
%
% Syntax:
%  [outputRay,rayPath, angleError] = calcNodalRay(eye,G,X,cameraMedium)
%
% Description
%   Given an eye structure and a coordinate on the retinal surface, the
%   routine returns a matrix that contains the path of a ray that satisfies
%   the property that the angle (wrt the optical axis) of the ray as it
%   departs the retina is equal to the angle with which it departs the
%   cornea. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   G                     - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X                     - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
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
    [outputRay,rayPath, angleError] = calcNodalRay(eye,[],eye.landmarks.vertex.coords)
%}


if nargin<2
    error('Invalid number of input arguments');
end

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
opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName','retinaToCamera', 'cameraMedium', cameraMedium );

% Define some options for the fmincon call in the loop
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% Define an error function that reflects the difference in angles from the
% initial ray and the output ray from the optical system
myError = @(p) calcOffsetFromParallel(opticalSystem,quadric.anglesToRay(X,p(1),p(2)));

% Supply an x0 guess which is the ray that connects the retinal point with
% the center of the lens.
[~, angle_p1p2, angle_p1p3] = quadric.angleRays( [0 0 0; 1 0 0]', quadric.normalizeRay([X'; eye.lens.center-X']') );
x0 = [angle_p1p2 angle_p1p3];

% Perform the search
[inputRayAngles, angleError] = fmincon(myError,x0,[],[],[],[],[-180,-180],[180,180],[],opts);

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
