function [G,X,angleError] = findRetinaFieldPoint( eye, degField, cameraMedium )
% Retinal point at a specified visual field location w.r.t. optical axis
%
% Syntax
%  fovea = human.landmarks.fovea( eye )
%
% Description:
%   Calculates the position on the retina of a point that has the specified
%   visual field position, where visual field position is defined with
%   respect to the optical axis of the eye.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   fovea                 - Structure with the subfields degField,
%                           geodetic, and coords
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    degField = [5.8 3.0 0];
    [G,X,angleError] = findRetinaFieldPoint( sceneGeometry.eye, degField);
    [outputRay,rayPath] = calcNodalRay(eye,G);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
%}

if nargin==2
    cameraMedium = 'air';
end

% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% Set some fmincon options we will be using below
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% To set an x0 guess, identify the retinal point that is at the specified
% angle w.r.t. to the optical axis and the center of the aperture stop.
R = rayFromAngles([-5.4 0 0]',-degField(2),-degField(1));
x0 = quadric.intersectRay(S,R,eye.retina.side,eye.retina.boundingBox);

% Convert the x0 zero guess into ellipsoidal geodetic coordinates
g0 = quadric.cartToEllipsoidalGeo( x0, S );

% Set the elev

% Define the optical axis ray
opticalAxis = [0 1; 0 0; 0 0];

% The objective is to match the angles between the output nodal ray and the
% optical axis to the specified degField angles
myObj = @(G) sqrt(sum((degField(1:2).*[-1 -1]-wrapAngleRays(calcNodalRay(eye,G,[],cameraMedium),opticalAxis)).^2));

% Set the bounds
lb = [-90 -180 0]';
ub = [90 180 0]';

% Perform the search
[G,angleError] = fmincon(myObj, g0, [], [], [], [], lb, ub, [], opts);

% Obtain the coords coordinates of the fovea
X = quadric.ellipsoidalGeoToCart(G,S)';

end

%% LOCAL FUNCTION

% Wrap 
function angles = wrapAngleRays(R0,R1)
[~, angles(1), angles(2) ] = quadric.angleRays(R0,R1);
end

% Converts angles relative to the optical axis to a unit vector ray.
function R = rayFromAngles(p,angle_p1p2,angle_p1p3)
u = [1; tan(angle_p1p2); tan(angle_p1p3)];
R = quadric.normalizeRay([p, u]);
end
