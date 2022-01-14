function [theta,eccentricity] = calcRetinaCartTo2DPol(eye,X,P)
% Returns 2D polar coordinate for a point on the retinal surface
%
% Syntax:
%  [theta,eccentricity] = calcRetinaCartTo2DPol(eye,X,P)
%
% Description
%   We define a retinal polar coordinate space with angle theta as:
%         0 - nasal retinal meridian
%        90 - superior retinal meridian
%       180 - temporal retinal meridian
%       270 - inferior retinal meridian
%   and eccentricity in geodesic mm. By default, the origin of this space
%   is the fovea. This accomplishes a 2D embedding of the 3D retinal
%   surface.
%
%   This function returns the 2D, polar coordinate for a 3D, Cartesian
%   coordinate on the retinal surface.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   X                     - 3x1 vector that is a location on the retinal
%                           surface.
%   P                     - 3x1 vector that provides a coordinate on the
%                           retinal surface that is the origin of the 2D
%                           polar space. If undefined, the foveal
%                           coordinate in the eye structure is used.
%
% Outputs:
%   theta                 - Scalar, in degrees, in the domain 0-360.
%   eccentricity          - Scalar, in mm.
%
% Examples:
%{
    % Test of invertability of the conversion
    eye = modelEyeParameters();
    initialTheta = 225; initialEccen = 1.75;
    X = calcRetina2DPolToCart(eye,initialTheta,initialEccen);
    [theta,eccentricity] = calcRetinaCartTo2DPol(eye,X);
    assert(abs(initialEccen-eccentricity)<1e-2);
    assert(abs(initialTheta-theta)<1e-2);
%}

arguments
    eye (1,1) {isstruct}
    X (3,1) {mustBeNumeric}
    P (3,1) {mustBeNumeric} = eye.landmarks.fovea.coords'
end

% Extract the quadric vector for the retinal surface
S = eye.retina.S;

% Check that X and P are on the retinal surface
myS = quadric.vecToFunc(S);
if myS(X(1),X(2),X(3))> 0.01
    error('calcRetina2DPolToCart:invalidPoint','Coordinate X must lie on retinal surface')
end
if myS(P(1),P(2),P(3))> 0.01
    error('calcRetina2DPolToCart:invalidPoint','Reference point P must lie on retinal surface')
end

% If P and X are extremely close, return zero for theta and eccentricity
if max(P-X) < eps
    theta = 0; eccentricity = 0;
    return
end

% Obtain the geodesic distance between X and P
eccentricity = quadric.geodesic(S,[X,P]);

% Define the vector between X and P
V = quadric.normalizeRay([P, X-P]);

% Obtain the surface normal at P
N = quadric.surfaceNormal(S,P);

% Define the rotation matrix that rotates the surface normal to the
% cardinal plane
GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
    norm(cross(A,B)) dot(A,B)  0;
    0              0           1];
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G*inv(Fi);
a = N(:,2); b = [1;0;0];
U=UU(FFi(a,b), GG(a,b));

% Apply the rotation to the vector V
V(:,2) = U*V(:,2);

% Obtain the thea angle of the vector V in the xy plane
theta = wrapTo360(atan2d(V(3,2),V(2,2)));

end
