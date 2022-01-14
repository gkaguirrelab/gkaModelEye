function X = calcRetina2DPolToCart(eye,theta,eccentricity,P)
% Returns the 3D coordinate of a retinal point defined in 2D polar space
%
% Syntax:
%  X = calcRetina2DPolToCart(eye,polarAngle,eccentricity,P)
%
% Description
%   We define a retinal polar coordinate space with angle theta as:
%         0 - nasal retinal meridian
%        90 - superior retinal meridian
%       180 - temporal retinal meridian
%       270 - inferior retinal meridian
%   and eccentricity in geodesic mm. By default, the origin of this space
%   is the fovea. This accomplishes a 2D embedding of the 3D retinal
%   surface. This function returns the 3D, Cartesian point on the retinal
%   surface corresponding to a 2D, polar retinal point.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   theta                 - Scalar, in degrees.
%   eccentricity          - Scalar, in mm.
%   P                     - 3x1 vector that provides a coordinate on the
%                           retinal surface that is the origin of the 2D
%                           polar space. If undefined, the foveal
%                           coordinate in the eye structure is used.
%
% Outputs:
%   X                     - 3x1 vector that is the location on the retinal
%                           surface that corresponds to the given polar
%                           coordinate.
%
% Examples:
%{
    % Basic operation for theta = 0, eccentricity = 1 mm.
    eye = modelEyeParameters();
    X = calcRetina2DPolToCart(eye,0,1);
%}
%{
    % Show points on the retina at 5 mm geodesic distance from the fovea
    eye = modelEyeParameters();
    plotOpticalSystem(eye);
    view([-50,10])
    cmap = colormap('hsv');
    for dd=1:3
        c(:,dd) = interp1(1:256,cmap(:,dd),linspace(1,255,360));
    end
    for theta = 0:15:345
        X = calcRetina2DPolToCart(eye,theta,5);
        plot3(X(1),X(2),X(3),'*','Color',c(theta+1,:));
    end
%}


arguments
    eye (1,1) {isstruct}
    theta (1,1) {mustBeNumeric}
    eccentricity (1,1) {mustBeNumeric}
    P (3,1) {mustBeNumeric} = eye.landmarks.fovea.coords';
end

% Extract the quadric vector for the retinal surface
S = eye.retina.S;

% Check that P is on the retinal surface
myS = quadric.vecToFunc(S);
if myS(P(1),P(2),P(3))> 0.01
    error('calcRetina2DPolToCart:invalidPoint','Reference point P must lie on retinal surface')
end

% If the eccentricity is zero, return P
if eccentricity < eps
    X = P;
    return
end

% Define a vector in the XY plane that is at angle theta
V = [P, [0;cosd(theta);sind(theta)]];

% Obtain the surface normal at P
N = quadric.surfaceNormal(S,P);

% Define the rotation matrix that rotates the normal of the plane of V to
% P. Rotation matrix solution from:
%   https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/897677#897677
GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
    norm(cross(A,B)) dot(A,B)  0;
    0              0           1];
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G*inv(Fi);
a = [1;0;0]; b = N(:,2);
U=UU(FFi(a,b), GG(a,b));

% Apply the rotation to the vector V
V(:,2) = U*V(:,2);

% Define an objective, which is the distance between the desired and
% measured geodesic
myObj = @(d) eccentricity - geoFromEuc(d,V,S);

% Find the distance along vector V needed to produce a geodesic distance
% equal to the desired eccentricity
d = fzero(myObj,eccentricity);

% Find the retinal point corresponding to this solution
G = V(:,1)+d*V(:,2);
[~,X] = quadric.distancePointEllipsoid( G, S );

end

% Local function to return the geodesic distance given the distance along
% the vector V, projected to the closest point on quadric S
function geoDistance = geoFromEuc(d,V,S)
G = V(:,1)+d*V(:,2);
[~,X] = quadric.distancePointEllipsoid( G, S );
geoDistance = quadric.geodesic(S,[V(:,1),X]);
end
