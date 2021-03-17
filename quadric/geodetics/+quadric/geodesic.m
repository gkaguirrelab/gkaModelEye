function [distance,geodesicPathCoords] = geodesic(S,P,pathResolution,surfaceTol)
% Find the geodesic distance between points on a tri-axial ellipsoid
%
% Syntax:
%  distance = quadric.geodesic(S,p,pathResolution)
%
% Description:
%   The "inverse" geodesic problem identifies the minimum distance between
%   two points on the tri-axial ellipsoidal surface. There have been many
%   treatments of this problem, which vary in their accuracy and robustness
%   to the special conditions that occur in the vicinity of the umbilical
%   points of the ellipsoid. I am unaware of any general, exact solution
%   that is able to solve the inverse problem for any arbitrary pair of
%   points. A particular difficulty in the present application is that the
%   fovea lies close to the pole of the elliopsoidal surface, making it
%   challenging to calculate the exact geodesic around this location.
%
%   This routine provides an approximation to the solution by identifying a
%   plane that intersects the ellipsoid, through the two points, and
%   provides the minimum arc length between the points along the ellipse.
%   This plane will be similar to the "normal section curve" which closely
%   approximates the geodesic on an ellipsoid with rotational symmetry. For
%   pairs of points, the result appears to be accurate to within ~ 1 /
%   5,000 of the value provided by the Panou 2013 boundary solution (see:
%   geodesicPanou). If a third, intermediate point is provided, then the
%   geodesic will be constrained to pass through that point as well.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   P                     - 3x2 or 3x3 matrix that specifies the points 
%                           through which the geodesic will pass. The
%                           values can be given in Cartesian or geodetic
%                           coordinates. For the latter, the values are
%                           beta, omega, and elevation in units of degrees.
%                           Beta is defined over the range -90:90, and
%                           omega over the range -180:180. Elevation has an
%                           obligatory value of zero as this solution is
%                           only defined on the surface. The routine will
%                           attempt to determine which type of coordinate
%                           set has been provided.
%   pathResolution        - Scalar. The number of points for the 
%                           geodeticPathCoords
%   surfaceTol            - Scalar. When interpreting if the values in p
%                           are Cartesian or geodetic coordinates, this is
%                           the tolerance within which a candidate
%                           Cartesian coordinate must be on the quadric
%                           surface.
%
% Outputs:
%   distance              - 1x(n-1) vector. Geodesic distance from the 
%                           first point in p to each subsequent point.
%   geodesicPathCoords    - 3xpathResolution matrix of locations along the
%                           geodesic
%
% Examples:
%{
    % Distance from the fovea to the optic disc
    eye = modelEyeParameters('sphericalAmetropia',0);
    S = eye.retina.S;
    P = [eye.landmarks.fovea.coords',eye.landmarks.opticDisc.coords'];
    [distance,geodesicPathCoords] = quadric.geodesic(S,P);
    [distancePanou,geodesicPanau] = quadric.geodesicPanou(S,P);
    assert( abs(distance-distancePanou)/distancePanou < 5e-4);
    outline = sprintf('Geodesic distance (ellipse approximation) from the fovea to the optic disc: %2.2f mm\n',distance);
    fprintf(outline);
    % Plot the geodesic
    boundingBox = [-30 30 -30 30 -30 30];
    figure
    quadric.plotImplicitSurface(S, boundingBox, 'k', 0.25, 'red');
    camlight
    hold on
    plot3(geodesicPathCoords(1,:),geodesicPathCoords(2,:),geodesicPathCoords(3,:),'.k');
    plot3(geodesicPanau(1,:),geodesicPanau(2,:),geodesicPanau(3,:),'-r');
    X0 = eye.landmarks.fovea.coords;
    X1 = eye.landmarks.opticDisc.coords;
    plot3(X0(1),X0(2),X0(3),'*r');
    plot3(X1(1),X1(2),X1(3),'*b');
%}


arguments
    S {mustBeNumeric}
    P (3,:) {mustBeNumeric}
    pathResolution (1,1) {mustBeNumeric} = 50
    surfaceTol (1,1) {mustBeNumeric} = 1e-6
end

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% p should be either size (3:2) or (3:3)
if size(P,2)<2 || size(P,3)>3
        error('geodesic:invalidCoordinate','Supply either 2 or 3 surface points in p.')
end

% Interpret p and convert to Cartesian coordinates if needed.
Sfunc = quadric.vecToFunc(S);
if Sfunc(P(1,1),P(2,1),P(3,1)) > surfaceTol
    
    % If the last value of the coordinate is not zero, then this can't be a
    % geodetic coordinate either
    if P(3,1) > surfaceTol
        error('geodesic:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the surface of S.')
    end
    
    % Check that the candidate beta and omega values are in range
    if any(abs(P(1,:))>90) || any(abs(P(2,:))>180)
        error('geodesic:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the surface of S.')
    end
    
    % Looks like valid geodetic coordinates. Convert to Cartesian
    for ii=1:size(P,2)
        P(:,ii) = quadric.ellipsoidalGeoToCart( P(:,ii), S );
    end
end

% Bounds
bound = max(quadric.radii(S));

% Objective
myObj = @(x) objective(x,P,S,pathResolution);

% Options
options = optimoptions('fmincon');
options.Display = 'off';

% Search
[x1,fval1]=fminbnd(myObj,-bound,bound);
[x2,fval2]=fmincon(myObj,0,[],[],[],[],-bound,bound,[],options);

if fval1<fval2
    x = x1;
else
    x = x2;
end

% Call again to obtain geodeticPathCoords
[distance,geodesicPathCoords] = objective(x,P,S,pathResolution);

end


%% LOCAL FUNCTIONS

function [d,geodeticPathCoords] = objective(x,P,S,pathResolution)
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
Pc = P-quadricCenter;

if size(Pc,2)==2
    % Define a plane with the reflection of the initial points.
    % Find the normal to this initial plane, and then adjust by x units
    P3 = -mean(Pc,2);
    u = cross(Pc(:,1)-Pc(:,2),Pc(:,1)-P3);
    Pc(:,3) = P3+x.*u;
end

% Parameters of the plane equation
planeNormal=cross(Pc(:,1)-Pc(:,2),Pc(:,1)-Pc(:,3));
A=planeNormal(1);B=planeNormal(2);C=planeNormal(3);
D = -dot(planeNormal,Pc(:,1));

% The intersection of the plane and the ellipsoid defines an ellipse. Aye
% Bye are the semi-major and semi-minor axis lengths, and [q1; q2; q3] is
% the ellipse center
r = quadric.radii(Sc);
a=r(1);b=r(2);c=r(3);
[Aye,Bye,q1,q2,q3] = quadric.intersectPlaneEllipsoid(A,B,C,D,a,b,c);
Ec = [q1;q2;q3];

% The ellipse lies in the plane identifed by Pc. For the calculations that
% follow, it is easier if the ellipse is parallel to the xy plane. To do
% so, we rotate the points PC, and the ellipse center, eC, so that these
% are all parallel to the xy plane.
w = planeNormal/norm(planeNormal);
R = [null(w'),w];
if det(R)<0, R(:,1:2) = R(:,2:-1:1); end
PcR = (Pc'*R)';
EcR = (Ec'*R)';

% p1 and p2 are the X0 and X1 coordinates, now with reference to the
% ellipse plane
PcRc = PcR - EcR;

% Find the angle of each of these points with respect to the ellipse
% center.
t1 = atan2(PcRc(2,1),PcRc(1,1));
t2 = atan2(PcRc(2,2),PcRc(1,2));

% Obtain the arc lengths around the ellipse to each of these points using
% elliptic integration of the second kind
k2=sqrt(1-Aye^2/Bye^2);
fun2=@(angle) sqrt(1-k2^2*(sin(angle)).^2);
arc1=Bye*integral(fun2,0,t1);
arc2=Bye*integral(fun2,0,t2);

% This is the objective.
d = abs(arc1-arc2);

% If requested, generate points on the geodesic
if nargout == 2
    thetas = linspace(t1,t2,pathResolution);
    tmp = mod(thetas+pi/2,pi);
    thetas = tmp+pi*(thetas>0&tmp==0)-pi/2;
    
    % Figure out the orientation of Aye and Bye with respect to the x and y
    % dimensions
    if abs(mean(PcRc(1,:).^2/(Bye^2) +  PcRc(2,:).^2/(Aye^2) - 1)) > 5e-2
        Xf = @(theta) (Aye.*Bye)./(Aye.^2*tan(theta).^2 + Bye.^2).^(1/2);
        Yf = @(theta) (Aye.*Bye.*tan(theta))./(Aye.^2*tan(theta).^2 + Bye.^2).^(1/2);
    else
        Xf = @(theta) (Aye.*Bye)./(Bye.^2*tan(theta).^2 + Aye.^2).^(1/2);
        Yf = @(theta) (Aye.*Bye.*tan(theta))./(Bye.^2*tan(theta).^2 + Aye.^2).^(1/2);
    end

    % Get the points from the thetas
    Xp = Xf(thetas);
    Yp = Yf(thetas);

    % Convert back to the space of S
    geodeticPathCoords = (([Xp; Yp; zeros(size(Xp))] + EcR)'*R')' + quadricCenter;    
    
end

end

