function [G,X,angleError] = findRetinaFieldPoint( eye, degField, cameraMedium )
% Retinal point at a specified visual field location w.r.t. optical axis
%
% Syntax
%  [G,X,angleError] = findRetinaFieldPoint( eye, degField, cameraMedium )
%
% Description:
%   Calculates the position on the retina of a point that has the specified
%   visual field position, where visual field position is defined with
%   respect to the optical axis of the eye. The identified retinal point
%   has the property that a nodal ray that emerges from this point exits
%   the cornea at angles with respect to the optical axis that are equal to
%   the values specified in degField.
%
% Inputs
%   eye                   - Structure.
%   degField              - 2x1 vector that specifies the visual field
%                           position of a point in degree angles along the
%                           horizontal and vertical meridian with respect
%                           to the optical axis of the eye.
%   cameraMedium          - Char vector. Used to determine the refractive
%                           index of the medium of the visual field.
%                           Defaults to "air" if not passed.
%
% Outputs
%   G                     - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X                     - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   angleError            - Scalar. The angle between the desired output
%                           ray and the solution. 
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    degField = [5.8 3.0];
    [G,X,angleError] = findRetinaFieldPoint( sceneGeometry.eye, degField);
    [outputRay,rayPath] = calcNodalRay(sceneGeometry.eye,G);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
%}

if nargin<2
    error('Need to specify an eye structure and the visual field angles');
end

if nargin==2
    cameraMedium = 'air';
end

% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% Set some fmincon options we will be using below
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% To set an x0 guess, identify the retinal point that is at the specified
% angle w.r.t. to the optical axis and the center of the aperture stop.
R = quadric.anglesToRay(eye.stop.center',degField(1),degField(2));
x0 = quadric.intersectRay(S,R,eye.retina.side,eye.retina.boundingBox);

% Convert the x0 zero guess into ellipsoidal geodetic coordinates, beta
% (latitude), omega (longitude), and elevation.
g0 = quadric.cartToEllipsoidalGeo( x0, S );

% Force elevation to be zero
g0(3) = 0;

% Define the optical axis ray
opticalAxis = [0 1; 0 0; 0 0];

% The objective is to match the angles between the output nodal ray and the
% optical axis to the specified degField angles.  Note that the elevational
% angle is inverted. This is because a negative value in this context
% corresponds to deflection of the visual axis upwards in the visual field.
myObj = @(G) sqrt(sum((-degField(1:2).*[1 -1]-wrapAngleRays(calcNodalRay(eye,G,[],cameraMedium),opticalAxis)).^2));

% Set the bounds
lb = [-90 -180 0]';
ub = [90 180 0]';

% Perform the search
[G,angleError] = fmincon(myObj, g0, [], [], [], [], lb, ub, [], opts);

% Obtain the Cartesian coordinates of the fovea
X = quadric.ellipsoidalGeoToCart(G,S);

end

%% LOCAL FUNCTION

% Wrap 
function angles = wrapAngleRays(R0,R1)
[~, angles(1), angles(2) ] = quadric.angleRays(R0,R1);
end
