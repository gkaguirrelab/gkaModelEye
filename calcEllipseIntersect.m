function [intersectionCoords, curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
% Returns coords and curvature of an ellipse intersected by a ray
%
% Syntax:
%  [intersectionCoords, curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
%
% Description:
%   Implements trigonometric operations to identify the point at which a
%   line intersects an ellipse, and returns the curvature of the ellipse at
%   the point of contact.
%
% Inputs:
%   coordsInitial         - 2x1 vector, with the values corresponding to
%                           the z-position and height of the initial
%                           position of the ray.
%   theta                 - Scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   ellipseCenterZ        - Scalar in mm. The position of the center of the
%                           elliptical surface on the z-axis
%   ellipseRadii          - A 2x1 vector in mm, with the values
%                           corresponding to the radii of the ellipitical
%                           surface along the z and height axes,
%                           respectively.
%
% Outputs:
%   intersectionCoords    - 2x1 vector with the values corresponding to the
%                           z-position and height at which the ray strikes
%                           the elliptical surface
%   curvature             - Scalar. The radius of curvature of the
%                           elliptical surface at the point of intersection
%                           with the ray.
%
% Examples:
%{
    [intersectCoords,curvature] = calcEllipseIntersect([-3.7,0], 0, -9.1412, [9.1412, 8.4277] )
%}


% Initialize return variables
intersectionCoords = [nan, nan];
curvature = nan;

% Store the sign of the radius values. They radii must have the same sign.
if prod(sign(ellipseRadii))==-1
    error('rayTraceCenteredSurfaces:incompatibleConvexity','The radii of the elliptical lens surface must have the same sign.');
end
radiiSign = sign(ellipseRadii(1));
% Convert the radii to their absolute values
ellipseRadii = abs(ellipseRadii);

% Convert the ray position and theta to the slope (m) and y-axis
% intercept(c) of a line.
m = tan(theta);
c = coordsInitial(2)-m*coordsInitial(1);

% Obtain the pair of z and h coordinates at which the line will intersect
% the ellipse
M = (1/ellipseRadii(1)^2) + (1/ellipseRadii(2)^2)* m^2;
N = 2*(1/ellipseRadii(2)^2)*m*c - 2*(1/ellipseRadii(1)^2)*ellipseCenterZ;
O = (1/ellipseRadii(2)^2)*c^2 + (1/ellipseRadii(1)^2)*ellipseCenterZ^2 -1 ;
determinant = (N^2 - 4* M * O);

% Detect if we have a tangential or non-intersecting ray
if ~isreal(determinant) || determinant<0
    return
end

zCoords = [(- N + sqrt(determinant))/ (2*M), (- N - sqrt(determinant))/ (2*M)];
hCoords = [m*zCoords(1) + c, m*zCoords(2) + c];

% If the radiiSign is positive, report the coordinates on the left-hand
% side of the ellipse, otherwise report the coordinates on the right
if radiiSign<0
    [~,rightIdx] = max(zCoords);
    intersectionCoords = [zCoords(rightIdx),hCoords(rightIdx)];
else
    [~,leftIdx] = min(zCoords);
    intersectionCoords = [zCoords(leftIdx),hCoords(leftIdx)];
end

% Detect if we have a tangential or non-intersecting ray
P = (intersectionCoords(1)-ellipseCenterZ)/ellipseRadii(1);
if ~isreal(P)
    return
end
if P>1 || P<-1
    return
end

% Calculate the radius of curvature at the point of intersection. If
% radiiSign is negative, report a negative radius of curvature
t = acos(P);
curvature = radiiSign*((ellipseRadii(1)^2*sin(t)^2 + ellipseRadii(2)^2*cos(t)^2)^(3/2))/(ellipseRadii(1)*ellipseRadii(2));

end