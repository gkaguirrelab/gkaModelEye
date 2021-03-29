function p = plotRay(R,t,lineColor,lineAlpha)
% Add a 3D plot of a ray to the active figure
%
% Syntax:
%  p = quadric.plotRay(R,t,lineColor,lineAlpha)
%
% Description:
%   Plot a ray, scaled by t, using the supplied line color and alpha.
%
% Inputs:
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%   t                     - Scalar. A multiplier of the u component of R.
%   lineColor             - 1x3 vector or char string that specifies a
%                           valid MATLAB plor color. E.g. [0.5 0.5 0.5] or
%                           'red'.
%   lineAlpha             - Scalar, range 0-1. Specifies the transparency
%                           of the line from 0 (opaque) to 1
%                           (fully transparent)
%
% Outputs:
%   p                     - Handle to the line object
%
% Examples:
%{
    
%}

% Handle incomplete input arguments
if nargin==1
    t = 1;
    lineColor=[0.9 0 00];
    lineAlpha=0.8;
end

if nargin==2
    lineColor=[0.9 0 00];
    lineAlpha=0.8;
end

if nargin==3
    lineAlpha=0.8;
end

% Get the points
xx = [R(1,1),R(1,1)+R(1,2)*t];
yy = [R(2,1),R(2,1)+R(2,2)*t];
zz = [R(3,1),R(3,1)+R(3,2)*t];

% Plot the line
p = plot3(xx,yy,zz,'-','Color',lineColor);
p.Color(4) = lineAlpha;


end