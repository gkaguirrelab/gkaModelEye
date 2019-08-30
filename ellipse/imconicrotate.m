function pr = imconicrotate(p, theta)
% Rotates a conic section given as an implicit equation.
% The conic section is assumed to be centered at the origin.
%
% Input arguments:
% p:
%    the parameter vector p = [a b c d f g]
% theta:
%    the angle to rotate
%
% See also: imconicrotation, imconictranslate

% Copyright 2010 Levente Hunyadi

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
assert(length(p) == 6, 'imconicrotate:DimensionMismatch', ...
    'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
validateattributes(theta, {'numeric'}, {'real','scalar'});

a = p(1);
b = 0.5 * p(2);
c = p(3);
d = 0.5 * p(4);
f = 0.5 * p(5);
g = p(6);

ar = a*cos(theta)^2 - 2*b*cos(theta)*sin(theta) + c*sin(theta)^2;
br = b*(cos(theta)^2 - sin(theta)^2) + (a-c)*sin(theta)*cos(theta);
cr = a*sin(theta)^2 + 2*b*sin(theta)*cos(theta) + c*cos(theta)^2;
dr = d*cos(theta) - f*sin(theta);
fr = d*sin(theta) + f*cos(theta);

pr = zeros(size(p));
pr(:) = [ar 2*br cr 2*dr 2*fr g];
