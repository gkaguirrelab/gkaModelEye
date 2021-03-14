function distance = geodesic(S,G0,G1,X0,X1)
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
%   points of the ellipsoid.
%
%   This routine provides an approximation to the solution by identifying
%   a plane that intersects the ellipsoid, passes through the two points,
%   and provides the minimum arc length between the points alon the
%   ellipse. The result is accurate to within ~ 1 / 10,000 of the value
%   provided by the Panou 2013 boundary solution.
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
%
% Outputs:
%   distance              - Scalar. Distance of the geodetic between the
%                           two points.
%
% Examples:
%{
    % Distance from the fovea to the optic disc
    eye = modelEyeParameters();
    S = eye.retina.S;
    G0 = eye.landmarks.fovea.geodetic;
    G1 = eye.landmarks.opticDisc.geodetic;
    panouDistance = quadric.panouGeodesicDistance(S,G0,G1);
    distance = quadric.geodesic(S,G0,G1);
    assert( abs(distance-panouDistance)/panouDistance < 1e-3);
    outline = sprintf('Geodetic distance from the fovea to the optic disc: %2.2f mm\n',distance);
    fprintf(outline);
%}


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Convert geodetic coordinates to Cartesian
if nargin==3
    % Obtain the ellipsoidal geodetic coordinates for the two points, and
    % convert to radians
    X0 = quadric.ellipsoidalGeoToCart( G0, S );
    X1 = quadric.ellipsoidalGeoToCart( G1, S );
end

% Bounds
bound = max(quadric.radii(S));

% Objective
myObj = @(x) objective(x,X0,X1,S);

% Search
[~,distance]=fminbnd(myObj,-bound,bound);


end

%% LOCAL FUNCTIONS

function d = objective(x,X0,X1,S)
% Returns the arc distance along an ellipse on the ellipsoid surface


% Center the quadric
c = quadric.center(S);
Sc = quadric.translate(S,-c);

% Adjust the points for the center translation
X0c = X0-c;
X1c = X1-c;

% Create a third point that is the reflection across the quadric of the
% midpoint between X1c and X2c
X2c = -(X0c+X1c)/2;

% Find the normal to this initial plane, and then adjust X2c by x units.
X2c = X2c+x.*cross(X0c-X1c,X0c-X2c);

% Parameters of the plane equation
cp=cross(X0c-X1c,X0c-X2c);
A=cp(1);B=cp(2);C=cp(3);
D = -dot(cp,X0c);

% Intersection of the plane and the ellipsoid. Aye Bye are the semi-major
% and semi-minor axis lengths, and [q1; q2; q3] is the ellipse center
r = quadric.radii(Sc);
a=r(1);b=r(2);c=r(3);
[Aye,Bye,q1,q2,q3]=EllipsoidPlaneIntersection(A,B,C,D,a,b,c);
eC = [q1;q2;q3];

% The ellipse lies in the plane identifed by X0c, X1c, X2c. For the
% calculations that follow, it is easier if the ellipse is parallel to the
% xy plane. To do so, we rotate the points X0c, X1c, X2c, and the ellipse
% center ec, so that these are all parallel to the xy plane.
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
t1 = atan2d(p1(2),p1(1));
t2 = atan2d(p2(2),p2(1));

% Obtain the arc lengths around the ellipse to each of these points using
% elliptic integration of the second kind
k2=sqrt(1-Aye^2/Bye^2);
fun2=@(angle) sqrt(1-k2^2*(sin(angle)).^2);
arc1=Bye*integral(fun2,0,deg2rad(t1));
arc2=Bye*integral(fun2,0,deg2rad(t2));

% This is the objective.
d = abs(arc1-arc2);

end


%% Functions by Sebahattin Bektas

function[Aye,Bye,q1,q2,q3]=EllipsoidPlaneIntersection(A,B,C,D,a,b,c)
% This function computes the intersection an Ellipsoid and a Plane,
% This source code was adapted from the implementation of the algorithms described in
%   P. Klein, "On the Ellipsoid and Plane Intersection Equation,"
%Applied Mathematics, Vol. 3 No. 11, 2012, pp. 1634-1640. doi: 10.4236/am.2012.311226.

% Plane equation          A.x + B.y + C.z + D = 0
% Ellipsoid equation    (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
%
% Author: Sebahattin Bektas,Ondokuz Mayis University,2015
%
% Inputs -----------------------------------------------------------------
%   A,B,C,D    :Plane equation's coefficients
%
%   a,b,c      : Semi-axes of ellipsoid
%
% Outputs ----------------------------------------------------------------
%   Aye:  Semi-major axis of intersection ellipse
%   Bye:  Semi-minor axis of intersection ellipse
%
% q1,q2,q3 : Cartesian coordinates of intersection ellipse's center
% This source code is given for free! However, I would be grateful if you refer
% to corresponding article in any academic publication that uses
% this code or part of it. Here are the corresponding references:
%  BEKTAS, Sebahattin. Orthogonal distance from an ellipsoid. Bol. Ciênc. Geod. [online]. 2014,
%  vol.20, n.4, pp. 970-983. ISSN 1982-2170.
%
car=A*B*C*D;
if car~=0
    kx2=1/a^2 + (A^2)/(C^2*c^2);
    ky2=1/b^2 + (B^2)/(C^2*c^2);
    kxy=(2*A*B)/(C^2*c^2);
    kx=+ (2*A*D)/(C^2*c^2);
    ky=+ (2*B*D)/(C^2*c^2) ;
    ksab=D^2/(C^2*c^2)- 1;
    ParA=[kx2  kxy  ky2  kx  ky  ksab];
    G=AtoG(ParA);
    q1=G(1);q2=G(2);q3=(A*q1+B*q2+D)/-C;
end
if C==0 & A~=0
    ParA=[(1/b^2 + (B/A/a)^2) 0  1/c^2 (2*D*B/(A*a)^2)  0  (-1+(D/A/a)^2)];G=AtoG(ParA);
    q2=G(1);q3=G(2);q1=(D+B*q2+C*q3)/-A;
end
if C==0 & B~=0
    ParA=[(1/a^2 + (A/B/b)^2) 0  1/c^2 (2*D*A/(B*b)^2)  0  (-1+(D/B/b)^2)];G=AtoG(ParA);
    q1=G(1);q3=G(2);q2=(D+A*q1+C*q3)/-B;
end
if A==0 & B==0 ,    q1=0;q2=0;q3=-D/C;end
if D==0,    q1=0;q2=0;q3=0;
end

quz=sqrt(q1^2+q2^2+q3^2);
n1=A/sqrt(A^2+B^2+C^2);
n2=B/sqrt(A^2+B^2+C^2);
n3=C/sqrt(A^2+B^2+C^2);
kap=(q1*A+B*q2+C*q3)/sqrt(A^2+B^2+C^2);
d= kap^2*(A^2+B^2+C^2)/(a^2*A^2+b^2*B^2+c^2*C^2);
ak=1;
bk=-(n1^2*(1/b^2+1/c^2)+n2^2*(1/c^2+1/a^2)+n3^2*(1/b^2+1/a^2));
ck=(n1/b/c)^2+(n2/a/c)^2+(n3/b/a)^2;
p=[ak bk ck];
kok=roots(p);
Bye=sqrt((1-d)/kok(1));
Aye=sqrt((1-d)/kok(2));
end


% AtoG
% Convert conic section equation from ABCDEF form to {center, axes, angle}.
% Minor modifications made to the original AtoG posted by Hui Ma
%  Note: Ma claimed copyright, but that is automatically released to the
% public upon his posting to FileExchange
%
%  ParA = [A,B,C,D,E,F]'-  parameter vector of the conic:  Ax^2 + Bxy + Cy^2 +Dx + Ey + F = 0
%  to geometric parameters  ParG = [Center(1:2), Axes(1:2), Angle]'
%
% The Angle value is in radians
%  code is:  1 - ellipse
%            2 - hyperbola
%            3 - parabola
%           -1 - degenerate cases
%            0 - imaginary ellipse
%            4 - imaginary parelle lines
%
%
function [ParG,code] = AtoG(ParA)
tolerance1 = 1.e-10;
tolerance2 = 1.e-20;
% ParA = ParA/norm(ParA);
if (abs(ParA(1)-ParA(3)) > tolerance1)
    Angle = atan(ParA(2)/(ParA(1)-ParA(3)))/2;
else
    Angle = pi/4;
end
c = cos(Angle);  s = sin(Angle);
Q = [c s; -s c];
M = [ParA(1)  ParA(2)/2;  ParA(2)/2  ParA(3)];
D = Q*M*Q';
N = Q*[ParA(4); ParA(5)];
O = ParA(6);
if (D(1,1) < 0) && (D(2,2) < 0)
    D = -D;
    N = -N;
    O = -O;
end
UVcenter = [-N(1,1)/2/D(1,1); -N(2,1)/2/D(2,2)];
free = O - UVcenter(1,1)*UVcenter(1,1)*D(1,1) - UVcenter(2,1)*UVcenter(2,1)*D(2,2);
% if the determinant of [A B/2 D/2;B/2 C E/2;D/2 E/2 F]is zero
% and if K>0,then it is an empty set;
% otherwise the conic is degenerate
Deg =[ParA(1),ParA(2)/2,ParA(4)/2;...
    ParA(2)/2,ParA(3),ParA(5)/2;...
    ParA(4)/2,ParA(5)/2,ParA(6)];
K1=[ParA(1),ParA(4)/2;ParA(4)/2 ParA(6)];
K2=[ParA(3),ParA(5)/2;ParA(5)/2 ParA(6)];
K = det(K1)+ det(K2);
if (abs(det(Deg)) < tolerance2)
    if (abs(det(M))<tolerance2) &&(K > tolerance2)
        code = 4;  % empty set(imaginary parellel lines)
    else
        code = -1; % degenerate cases
    end
else
    if (D(1,1)*D(2,2) > tolerance1)
        if (free < 0)
            code = 1; % ellipse
        else
            code = 0; % empty set(imaginary ellipse)
        end
    elseif (D(1,1)*D(2,2) < - tolerance1)
        code = 2;  % hyperbola
    else
        code = 3;  % parabola
    end
end
XYcenter = Q'*UVcenter;
Axes = [sqrt(abs(free/D(1,1))); sqrt(abs(free/D(2,2)))];
if code == 1 && Axes(1)<Axes(2)
    AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
    Angle = Angle + pi/2;
end
if code == 2 && free*D(1,1)>0
    AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
    Angle = Angle + pi/2;
end
% some people never learn...
Angle = mod(Angle,pi);
% while Angle > pi
%     Angle = Angle - pi;
% end
% while Angle < 0
%     Angle = Angle + pi;
% end
ParG = [XYcenter; Axes; Angle];
end

