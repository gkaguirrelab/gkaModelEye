function [ Xp, Yp ] = ellipsePerimeterPoints( transparentEllipseParams, steps, phase, noise, method )
% Returns a set of points on the boundary of a transparent ellipse
%
% Syntax:
%  [ Xp, Yp ] = ellipsePerimeterPoints( transparentEllipseParams, steps )
%
% Description:
%   Given an ellipse in transparent form, the routine returns the
%   coordinates [Xp, Yp], which are points on the ellipse. The points are
%   distributed evenly by angle with respect to the center of the ellipse.
%   The number of returned points is specified in the input variable
%   "steps". Note that 5 points are needed to uniquely specify an ellipse,
%   and 6 are needed to characterize any error in fitting.
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
%   noise                 - Scalar. The width of the Gaussian distributed
%                           noise (in pixels) to be added to the X and Y
%                           positions of the perimeter points.
%   method                - Char vector. The method to be used to space the
%                           points about the ellipse perimeter. The choice
%                           renders the points equidistant in terms of:
%                             - 'arc': The elliptical arc length
%                             - 'angle': The angle with respect to the
%                                       center of the ellipse.
%
% Outputs:
%   Xp, Yp                - Each is a stepsx1 vector, providing the X
%                           and Y coordinate of each point.
%
% Examples:
%{
	transparentEllipseParams = [ 278   223   1546   0.23   1.9 ];
    [ Xp, Yp ] = ellipsePerimeterPoints(transparentEllipseParams,100,0,0,'arc');
    plot(Xp,Yp,'.r')
%}


% If steps was not passed, set to 5.
if nargin == 1
    steps = 5;
    phase = 0;
    noise = 0;
    method = 'angle';
end

if nargin == 2
    phase = 0;
    noise = 0;
    method = 'angle';
end

if nargin == 3
    noise = 0;
    method = 'angle';
end

if nargin == 4
    method = 'angle';
end

% Convert the transparent ellipse to explicit form
p = ellipse_transparent2ex(transparentEllipseParams);

% Distribute the parameters of the ellipse to individual variables
x = p(1);
y = p(2);
a = p(3);
b = p(4);
theta = p(5);

% Define the angles of the points with respect to the ellipse center.
switch method
    case 'angle'
        alpha = linspace(phase, 2*pi-(2*pi/steps)+phase, steps)';
    case 'arc'
        k1=sqrt(1-b^2/a^2);
        fun1=@(angle) sqrt(1-k1^2*(sin(angle)).^2);
        d=integral(fun1,0,2*pi)/steps;
        alpha = zeros(steps,1);
        for pp = 1:steps-1
            fun2 = @(angle) integral(fun1,alpha(pp),angle)-d;
            alpha(pp+1) = fzero(fun2,alpha(pp));
        end
    otherwise
        error('Not a recognized method for spacing ellipse perimeter points');
end
        
% Perform the calculation
sintheta = sin(theta);
costheta = cos(theta);
sinalpha = sin(alpha);
cosalpha = cos(alpha);
Xp = x + (a * cosalpha * costheta - b * sinalpha * sintheta);
Yp = y + (a * cosalpha * sintheta + b * sinalpha * costheta);

% Add noise
if noise > 0
    Xp = Xp + randn(size(Xp)).*noise;
    Yp = Yp + randn(size(Yp)).*noise;
end

end % function -- ellipsePerimeterPoints

