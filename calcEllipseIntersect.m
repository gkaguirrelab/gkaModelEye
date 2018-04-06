function [intersectionCoords,curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
% Returns coords and curvature of an ellipse intersected by a ray 
%
% Syntax:
%  [intersectCoords,curvature] = calcEllipseIntersect(initialCoords, theta, ellipseRadii )
%
% Description:
%
%
% Inputs:
%   coordsInitial         - A 2x1 vector, with the values corresponding to
%                           the z-position and height of the initial
%                           position of the ray.
%   theta                 - A scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   ellipseCenter         - Scaler. The position of the center of the
%                           elliptical surface on the z-axis
%   ellipseRadii          - A 2x1 vector, with the values corresponding to
%                           the radii of the ellipitical surface along the
%                           z and height axes, respectively.
%
% Outputs:
%
%
%{
    [intersectCoords,curvature] = calcEllipseIntersect([-3.7,0], 0, -9.1412, [9.1412, 8.4277] )
%}

% Check and interpret the inputs
radiiSign = sign(ellipseRadii(1));
ellipseRadii = abs(ellipseRadii);

% slope (m) and y-axis intercept(c) of a line, derived from initial
% position and the theta of the ray.
m = tan(theta);
c = coordsInitial(2)-m*coordsInitial(1);

% Parameter for ellipse

M = (1/ellipseRadii(1)^2) + (1/ellipseRadii(2)^2)* m^2;
N = 2*(1/ellipseRadii(2)^2)*m*c - 2*(1/ellipseRadii(1)^2)*ellipseCenterZ;
O = (1/ellipseRadii(2)^2)*c^2 + (1/ellipseRadii(1)^2)*ellipseCenterZ^2 -1 ;

determinant = (N^2 - 4* M * O);

xCoords = [(- N + sqrt(determinant))/ (2*M), (- N - sqrt(determinant))/ (2*M)];
yCoords = [m*xCoords(1) + c, m*xCoords(2) + c];


% If the radiiSign is positive, report the coordinates on the left-hand
% side of the ellipse, otherwise report the coordinates on the right
intersectionCoords = [xCoords((radiiSign/2)+1.5),yCoords((radiiSign/2)+1.5)];

t = acos((intersectionCoords(1)-ellipseCenterZ)/ellipseRadii(1));

% If the radiiSign is negative, report a negative radius of curvature
curvature = radiiSign*((ellipseRadii(1)^2*sin(t)^2 + ellipseRadii(2)^2*cos(t)^2)^(3/2))/(ellipseRadii(1)*ellipseRadii(2));

end