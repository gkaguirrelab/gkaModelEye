function [distance,geodesicPathCoords] = geodesic(S,P,pathResolution,surfaceTol)
% Geodesic distance between points on a tri-axial ellipsoid
%
% Syntax:
%  [distance,geodesicPathCoords] = geodesic(S,P,pathResolution,surfaceTol)
%
% Description:
%   The "inverse" geodesic problem identifies the minimum distance between
%   two points on the tri-axial ellipsoidal surface. There have been many
%   treatments of this problem, which vary in their accuracy and robustness
%   to the special conditions that occur in the vicinity of the umbilical
%   points of the ellipsoid. I am unaware of a general implementation for
%   the exact solution for the inverse problem for an arbitrary pair of
%   points. A particular difficulty in the present application is that the
%   fovea lies close to the pole of the elliopsoidal surface, making it
%   challenging to calculate the exact geodesic around this location.
%
%   This routine provides an approximation to the solution by identifying a
%   plane that intersects the ellipsoid, through the two points, and
%   provides the minimum arc length between the points along the ellipse.
%   This plane will be similar to the "normal section curve" which closely
%   approximates the geodesic on an ellipsoid with rotational symmetry. For
%   pairs of points, the result appears to be accurate to within 0.5% of
%   the value provided by the Panou 2013 boundary solution (see:
%   quadric.geodesicPanou and v001_geodesics.mlx).
%
%   If a third point on the ellipsoidal surface is provided, then the
%   routine will return the arc length between the first two points for an
%   ellipse that lies on a plane that minimizes the distance of the points
%   from the plane. This functionality may be used to construct a set of
%   points that are constrained to be on the same arc.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   P                     - 3x2 or 3x3 matrix that specifies the points 
%                           through which the geodesic will pass. The
%                           values can be given in Cartesian or geodesic
%                           coordinates. For the latter, the values are
%                           beta, omega, and elevation in units of degrees.
%                           Beta is defined over the range -90:90, and
%                           omega over the range -180:180. Elevation has an
%                           obligatory value of zero as this solution is
%                           only defined on the surface. The routine will
%                           attempt to determine which type of coordinate
%                           set has been provided.
%   pathResolution        - Scalar. The number of points for the 
%                           geodesicPathCoords
%   surfaceTol            - Scalar. When interpreting if the values in P
%                           are Cartesian or geodesic, this is the
%                           tolerance within which a candidate Cartesian
%                           coordinate must be on the quadric surface.
%
% Outputs:
%   distance              - Scalar. Approximation to the geodesic distance   
%                           between the first two points in P.
%   geodesicPathCoords    - 3xpathResolution matrix of locations along the
%                           geodesic.
%
% Examples:
%{
    % Distance from the fovea to the optic disc
    eye = modelEyeParameters('sphericalAmetropia',0);
    S = eye.retina.S;
    P = [eye.landmarks.fovea.coords',eye.landmarks.opticDisc.coords'];
    [distance,geodesicPathCoords] = quadric.geodesic(S,P);
    outline = sprintf('Geodesic distance (ellipse approximation) from the fovea to the optic disc: %2.2f mm\n',distance);
    fprintf(outline);
    % Plot the geodesic
    boundingBox = [-30 30 -30 30 -30 30];
    figure
    quadric.plotImplicitSurface(S, boundingBox, 'k', 0.25, 'red');
    camlight
    hold on
    plot3(geodesicPathCoords(1,:),geodesicPathCoords(2,:),geodesicPathCoords(3,:),'.k');
    plot3(P(1,1),P(2,1),P(3,1),'*r');
    plot3(P(1,2),P(2,2),P(3,2),'*b');
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

% Test here that the passed quadric is an ellipsoid
if ~strcmp(quadric.classify( S ),'ellipsoid')
    error('quadric:geodesic','Can only find the geodesic on an ellipsoid');
end

% P should be either size (3:2) or (3:3)
if size(P,2)<2 || size(P,3)>3
        error('geodesic:invalidCoordinate','Supply either 2 or 3 surface points in p.')
end

% Interpret P and convert to Cartesian coordinates if needed.
Sfunc = quadric.vecToFunc(S);
if Sfunc(P(1,1),P(2,1),P(3,1)) > surfaceTol
    
    % If the last value of the coordinate is not zero, then this can't be a
    % geodesic coordinate either
    if P(3,1) > surfaceTol
        error('geodesic:invalidCoordinate','Supply a Cartesian or geodesic coordinate that is on the surface of S.')
    end
    
    % Check that the candidate beta and omega values are in range
    if any(abs(P(1,:))>90) || any(abs(P(2,:))>180)
        error('geodesic:invalidCoordinate','Supply a Cartesian or geodesic coordinate that is on the surface of S.')
    end
    
    % Looks like valid geodesic coordinates. Convert to Cartesian
    for ii=1:size(P,2)
        P(:,ii) = quadric.ellipsoidalGeoToCart( P(:,ii), S );
    end
end

% Objective
myObj = @(x) objective(x,P,S,pathResolution);

% Search
[x,distance] = fminsearch(myObj,0);

% For very short distances, the non-linear result can be shorter than the
% Euclidean distance between the points due to numerical error. If we
% encounter this situation, set x to 0, which corresponds to the great
% ellipse (the path defined by a plane passing through the points and the
% center of the ellipsoid). This value will be very close to the geodesic.
if distance < norm(P(:,1)-P(:,2))
    x=0;
end

% Call again to obtain geodesicPathCoords
[distance,geodesicPathCoords] = objective(x,P,S,pathResolution);

end


%% LOCAL FUNCTIONS

function [d,geodesicPathCoords] = objective(x,P,S,pathResolution)
% Returns the arc distance along an ellipse on the ellipsoid surface
%
% Description:
%   P will specify two points on the ellipsoidal surface. We define a plane
%   by reflecting the midpoint of these points across the center of the
%   ellipsoid, and then translating the reflected point along the normal of
%   the initial plane by x units. The intersection of the plane with the
%   ellipsoid is obtained, and the arc length along the intersection
%   ellipse between P1 and P2 is calculated and returned. This routine
%   allows us to search over values of x to minimize the returned value of
%   d, and thus identify the approximate geodesic.
%
%   If a third point is provided in P, then this third point is used to
%   define the plane, and x is ignored.
%

% Center the quadric
quadricCenter = quadric.center(S);
Sc = quadric.translate(S,-quadricCenter);

% Adjust the points for the center translation
Pc = P-quadricCenter;

% If we have not been provided with a third point, use the passed x value
% to define a point and thus a plane.
if size(Pc,2)==2
    % Define a plane with the reflection of the initial points.
    % Find the normal to this initial plane, and then adjust by x units
    P3 = -mean(Pc,2);
    u = cross(Pc(:,1)-Pc(:,2),Pc(:,1)-P3);
    u = u/norm(u);
    Pc(:,3) = P3+x.*u;
end

% Parameters of the plane equation
planeNormal=cross(Pc(:,1)-Pc(:,2),Pc(:,1)-Pc(:,3));
planeNormal = planeNormal/norm(planeNormal);
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
% follow, we rotate to plane to place it (and the ellipse) parallel to one
% of the cardinal coordinate axes. To do so, we rotate the points Pc, and
% the ellipse center, Ec, so that these are all parallel to the xy, xz, or
% yz plane.

% Anonymous functions that, given vectors a and b, returns a rotation
% matrix which can be used to rotate a to be parallel to b.
ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
RU = @(a,b) eye(3) + ssc(cross(a,b)) + ...
     ssc(cross(a,b))^2*(1-dot(a,b))/(norm(cross(a,b))^2);
 
% Define a vector that is the cardinal axis that is closest to being
% aligned with the planeNormal, and match the sign of the direction of this
% axis vector with the planeNormal.
axisVector = [0; 0; 0];
axisIdx = find(abs(planeNormal)==max(abs(planeNormal)),1);
axisVector(axisIdx) = 1*sign(planeNormal(axisIdx));

% Obtain the rotation matrix that would make the planeNormal parallel to
% the axisVector.
R = RU(planeNormal,axisVector);
 
% Apply the rotation matrix to Pc and Ec.
PcR = R*Pc;
EcR = R*Ec;

% Adjust for the ellipse center
PcRc = PcR - EcR;

% Find the angle of each of these points with respect to the ellipse
% center. There is some machinery here to figure out which dimension may be
% dropped.
tanIdx = [1 2 3];
tanIdx = tanIdx(tanIdx ~= axisIdx);
t1 = atan2(PcRc(tanIdx(2),1),PcRc(tanIdx(1),1));
t2 = atan2(PcRc(tanIdx(2),2),PcRc(tanIdx(1),2));

% If both thetas are closer to pi than to zero, flip them
piFlipFlag = false;
if all(abs([t1 t2])>pi/2)
    t1 = t1+(-sign(t1)*pi);
    t2 = t2+(-sign(t2)*pi);
    piFlipFlag = true;
end

% Obtain the arc lengths around the ellipse to each of these points using
% elliptic integration of the second kind
k2=sqrt(1-Aye^2/Bye^2);
fun2=@(angle) sqrt(1-k2^2*(sin(angle)).^2);
arc1=Bye*integral(fun2,0,t1);
arc2=Bye*integral(fun2,0,t2);

% This is arc length between P1 and P2, and thus the the objective.
d = abs(arc1-arc2);

% If requested, generate points on the geodesic
if nargout == 2
    
    % Interpolate thetas between t1 and t2
    thetas = linspace(t1,t2,pathResolution);
    
    % Flip them back to the other side of the ellipse if needed
    if piFlipFlag
        thetas = thetas-sign(thetas)*pi;
    end
    
    % Figure out the orientation of Aye and Bye with respect to the x and y
    % dimensions
    if abs(mean(PcRc(tanIdx(1),:).^2/(Bye^2) +  PcRc(tanIdx(2),:).^2/(Aye^2) - 1)) > 5e-2
        Xf = @(theta) -(Aye.*Bye)./(Aye.^2*tan(theta).^2 + Bye.^2).^(1/2);
        Yf = @(theta) -(Aye.*Bye.*tan(theta))./(Aye.^2*tan(theta).^2 + Bye.^2).^(1/2);
    else
        Xf = @(theta) -(Aye.*Bye)./(Bye.^2*tan(theta).^2 + Aye.^2).^(1/2);
        Yf = @(theta) -(Aye.*Bye.*tan(theta))./(Bye.^2*tan(theta).^2 + Aye.^2).^(1/2);
    end

    % Assemble the coords from the thetas
    coords = nan(3,pathResolution);
    coords(tanIdx(1),:) = Xf(thetas);
    coords(tanIdx(2),:) = Yf(thetas);
    coords(axisIdx,:) = zeros(1,pathResolution);
    
    % Convert back to the space of S
    geodesicPathCoords =  R'*(coords + EcR) + quadricCenter;    
    
end

end

