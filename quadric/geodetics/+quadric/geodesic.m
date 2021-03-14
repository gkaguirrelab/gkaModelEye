function [distance,geodeticPathCoords] = geodesic(S,G0,G1,X0,X1,pathResolution)
% Find the geodesic distance between two points on a tri-axial ellipsoid
%
% Syntax:
%  distance = quadric.geodesic(S,G0,G1,X0,X1)
%
% Description:
%   The "inverse" geodesic problem identifies the minimum distance between
%   two points on the tri-axial ellipsoidal surface. There have been many
%   treatments of this problem, which vary in their accuracy and robustness
%   to the special conditions that occur in the vicinity of the umbilical
%   points of the ellipsoid. I am unaware of any general, exact solution
%   that is able to solve the inverse problem for any arbitrary pairs of
%   points. A particular difficulty in the present application is that the
%   fovea lies close to the pole of the elliopsoidal surface, making many
%   it challenging to calculate the exact geodesic around this location.
%
%   This routine provides an approximation to the solution by identifying a
%   plane that intersects the ellipsoid, passes through the two points, and
%   provides the minimum arc length between the points along the ellipse.
%   The result appears to be accurate to within ~ 1 / 5,000 of the value
%   provided by the Panou 2013 boundary solution (see: geodesicPanou).
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   G0, G1                - 3x1 vectors that provide the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%   pathResolution        - Scalar. The number of points for the 
%                           geodeticPathCoords
%
% Outputs:
%   distance              - Scalar. Distance of the geodetic between the
%                           two points.
%   geodeticPathCoords    - 3xpathResolution matrix of locations along the
%                           geodetic
%
% Examples:
%{
    % Distance from the fovea to the optic disc
    eye = modelEyeParameters('sphericalAmetropia',0);
    S = eye.retina.S;
    G0 = eye.landmarks.fovea.geodetic';
    G1 = eye.landmarks.opticDisc.geodetic';
    [distance,geodeticPathCoords] = quadric.geodesic(S,G0,G1);
    panouDistance = quadric.panouGeodesicDistance(S,G0,G1);
    assert( abs(distance-panouDistance)/panouDistance < 1e-3);
    outline = sprintf('Geodetic distance (ellipse approximation) from the fovea to the optic disc: %2.2f mm\n',distance);
    fprintf(outline);
    % Plot the geodesic
    boundingBox = [-30 30 -30 30 -30 30];
    figure
    quadric.plotImplicitSurface(S, boundingBox, 'k', 0.25, 'red');
    camlight
    hold on
    plot3(geodeticPathCoords(1,:),geodeticPathCoords(2,:),geodeticPathCoords(3,:),'.k');
    X0 = eye.landmarks.fovea.coords;
    X1 = eye.landmarks.opticDisc.coords;
    plot3(X0(1),X0(2),X0(3),'*r');
    plot3(X1(1),X1(2),X1(3),'*b');
%}


arguments
    S {mustBeNumeric}
    G0 = []
    G1 = []
    X0 (3,1) {mustBeNumeric} = quadric.ellipsoidalGeoToCart( G0, S )
    X1 (3,1) {mustBeNumeric} = quadric.ellipsoidalGeoToCart( G1, S )
    pathResolution (1,1) {mustBeNumeric} = 50
end

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Bounds
bound = max(quadric.radii(S));

% Objective
myObj = @(x) objective(x,X0,X1,S,pathResolution);

% Search
x=fminbnd(myObj,-bound,bound);

% Call again to obtain geodeticPathCoords
[distance,geodeticPathCoords] = objective(x,X0,X1,S,pathResolution);

end


%% LOCAL FUNCTIONS

function [d,geodeticPathCoords] = objective(x,X0,X1,S,pathResolution)
% Returns the arc distance along an ellipse on the ellipsoid surface
%
% Description:
%   We have two points (X0, X1). We define a plane by reflecting the
%   midpoint of X0, X1 across the center of the ellipsoid, and then
%   translating this point along the normal of this initial plane by x
%   units. The intersection of this plane with the ellipsoid is obtained,
%   and the arc length along the intersection ellipse between X0 and X1 is
%   calculated and returned. This routine allows us to search over values
%   of x that minimize the returned value of d, and thus identify the
%   approximate geodesic.

% Center the quadric
quadricCenter = quadric.center(S);
Sc = quadric.translate(S,-quadricCenter);

% Adjust the points for the center translation
X0c = X0-quadricCenter;
X1c = X1-quadricCenter;

% Create a third point that is the reflection across the quadric of the
% midpoint between X1c and X2c
X2c = -(X0c+X1c)/2;

% Find the normal to this initial plane, and then adjust X2c by x units.
X2c = X2c+x.*cross(X0c-X1c,X0c-X2c);

% Parameters of the plane equation
cp=cross(X0c-X1c,X0c-X2c);
A=cp(1);B=cp(2);C=cp(3);
D = -dot(cp,X0c);

% The intersection of the plane and the ellipsoid defines an ellipse. Aye
% Bye are the semi-major and semi-minor axis lengths, and [q1; q2; q3] is
% the ellipse center
r = quadric.radii(Sc);
a=r(1);b=r(2);c=r(3);
[Aye,Bye,q1,q2,q3] = quadric.intersectPlaneEllipsoid(A,B,C,D,a,b,c);
eC = [q1;q2;q3];

% The ellipse lies in the plane identifed by X0c, X1c, X2c. For the
% calculations that follow, it is easier if the ellipse is parallel to the
% xy plane. To do so, we rotate the points X0c, X1c, X2c, and the ellipse
% center, eC, so that these are all parallel to the xy plane.
xyz = [X0c,X1c,X2c]';
w = cross(xyz(2,:)-xyz(1,:),xyz(3,:)-xyz(1,:));
w = w/norm(w);
R = [null(w),w.'];
if det(R)<0, R(:,1:2) = R(:,2:-1:1); end
xyz2 = xyz*R;
ec2 = eC'*R;

% p1 and p2 are the X0 and X1 coordinates, now with reference to the
% ellipse plane
p1 = xyz2(1,:)' - ec2';
p2 = xyz2(2,:)' - ec2';

% Find the angle of each of these points with respect to the ellipse
% center.
t1 = atan2(p1(2),p1(1));
t2 = atan2(p2(2),p2(1));

% Obtain the arc lengths around the ellipse to each of these points using
% elliptic integration of the second kind
k2=sqrt(1-Aye^2/Bye^2);
fun2=@(angle) sqrt(1-k2^2*(sin(angle)).^2);
arc1=Bye*integral(fun2,0,t1);
arc2=Bye*integral(fun2,0,t2);

% This is the objective.
d = abs(arc1-arc2);

% Generate some points that define the geodesic
thetas = linspace(t1,t2,pathResolution);
geodeticPathCoords = (([Bye.*cos(thetas); Bye.*sin(thetas); zeros(1,pathResolution)]' + ec2)*inv(R))'+quadricCenter;

end

