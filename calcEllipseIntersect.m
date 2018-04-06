function [intersectCoords,curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
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


%% FIXED AT ZERO
ellipseCenterH= 0; % center y? OBLIGATORY 0

% slope (m) and y-axis intercept(c) of a line, derived from initial
% position and the theta of the ray.
m = tan(theta);
c = coordsInitial(2)-m*coordsInitial(1);

% Parameter for ellipse

M = (1/ellipseRadii(1)^2) + (1/ellipseRadii(2)^2)* m^2;
N = 2*(1/ellipseRadii(2)^2)*m*c - 2*(1/ellipseRadii(1)^2)*ellipseCenterZ;
O = (1/ellipseRadii(2)^2)*c^2 + (1/ellipseRadii(1)^2)*ellipseCenterZ^2 -1 ;

determinant = (N^2 - 4* M * O);

x1  = (- N + sqrt(determinant))/ (2*M);
x2  = (- N - sqrt(determinant))/ (2*M);

y1 = m*x1 + c;
y2 = m*x2 + c;

cond1 = ( (     ((x1-ellipseCenterZ))^2 /ellipseRadii(1)^2  ) >= (1-1e-15)) && ...
    ( (((x1-ellipseCenterZ))^2/ellipseRadii(1)^2 ) <= (1+1e-15));

cond2 =  ( (((x2-ellipseCenterZ))^2/ellipseRadii(1)^2 ) >= (1-1e-15)) && ...
    ( (((x2-ellipseCenterZ))^2/ellipseRadii(1)^2  ) <= (1+1e-15));

if (cond1 == 1 && cond2 == 0)
    x2 = x1;
    y2 = y1;
elseif (cond1 == 0 && cond2 == 1)
    x1 = x2;
    y1 = y2;
elseif (cond1 == 0 && cond2 == 0)
    x1 = NaN;
    x2 = NaN;
    y1 = NaN;
    y2 = NaN;
end

% Need some logic here to decide if we return x1,y1 or x2,y2
intersectCoords = [x1,y1];

t = acos((x1-ellipseCenterZ)/ellipseRadii(1));

% radius of curvature
curvature = ((ellipseRadii(1)^2*sin(t)^2 + ellipseRadii(2)^2*cos(t)^2)^(3/2))/(ellipseRadii(1)*ellipseRadii(2));

end