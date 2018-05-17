function [ Xp, Yp ] = ellipsePerimeterPoints( transparentEllipseParams, steps, phase )
% Returns a set of points on the boundary of a transparent ellipse
%
% Syntax:
%  [ Xp, Yp ] = ellipsePerimeterPoints( transparentEllipseParams, steps )
%
% Description:
%   Given an ellipse in transparent form, the routine returns the
%   coordinates [Xp, Yp], which are points on the ellipse. The points are
%   distributed evenly by angle with respect to the center of the ellisep.
%   The number of returned points is specified in the input variable
%   "steps". Note that 5 points are needed to uniquely specify an ellipse.
%
% Inputs:
%   transparentEllipseParams - A 1x5 vector containing the parameters of an
%                           ellipse in transparent format.
%   steps                 - Scalar. The number of points on the boundary
%                          	to be returned. If not passed, steps is set to
%                          	5.
%   phase                 - Scalar. The phase (in radians) of the points
%                           around the pupil perimeter. If not passed,
%                           phase is set to zero.
%
% Outputs:
%   Xp, Yp                - Each is a stepsx1 vector, providing the X
%                           and Y coordinate of each point.
%

% If steps was not passed, set to 5.
if nargin == 1
    steps = 5;
    phase = 0;
end

if nargin == 2
    phase = 0;
end

% Convert the transparent ellipse to explicit form
p = ellipse_transparent2ex(transparentEllipseParams);

% Distribute the parameters of the ellipse to individual variables
x = p(1);
y = p(2);
a = p(3);
b = p(4);
theta = p(5);

% Perform the calculation
sintheta = sin(theta);
costheta = cos(theta);
alpha = linspace(phase, 2*pi-(2*pi/steps)+phase, steps)';
sinalpha = sin(alpha);
cosalpha = cos(alpha);
Xp = x + (a * cosalpha * costheta - b * sinalpha * sintheta);
Yp = y + (a * cosalpha * sintheta + b * sinalpha * costheta);

end % function -- ellipsePerimeterPoints

