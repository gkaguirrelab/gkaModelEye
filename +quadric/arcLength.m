function [d1, d2] = arcLength(S,p1,p2)
% Find the arc distance between two points on a quadric surface
%
% Syntax:
%  d = quadric.arcLength(S,p1,p2)
%
% Description:
%   Returns the arc length distance on the quadric surface between two
%   points. This currently only works for ellipsoids.
%
% Inputs:
%   S                     - 4x4 quadric surface matrix
%   p1, p2                - 3x1 vectors that specify the location of points
%                           on the quadric surface.
%
% Outputs:
%   d1, d2                - Scalars. Arc length distance.
%
% Examples:
%{
%}


% Pre-allocate the output variables
d = nan;

% Obtain the theta1 and theta2 values that define the angle of each point
% with respect to the center of the ellipsoid

ellipticIntegral_p1p2=@(theta) sqrt(1-sqrt(1-eye.posteriorChamber.radii(2).^2/eye.posteriorChamber.radii(1).^2)^2.*(sin(theta)).^2);
ellipticIntegral_p1p3=@(theta) sqrt(1-sqrt(1-eye.posteriorChamber.radii(3).^2/eye.posteriorChamber.radii(1).^2)^2.*(sin(theta)).^2);
arcLength_p1p2 = @(theta1,theta2) eye.posteriorChamber.radii(1).*integral(ellipticIntegral_p1p2, theta1, theta2);
arcLength_p1p3 = @(theta1,theta2) eye.posteriorChamber.radii(1).*integral(ellipticIntegral_p1p3, theta1, theta2);

% For the calculation, the first theta value is zero, as we are
% calculating distance from the posterior chamber apex (i.e., the
% intersection of the optical axis with the retina).
axes.visual.mmRetina = [arcLength_p1p2(0,deg2rad(axes.visual.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.visual.degRetina(2))), 0];
axes.opticDisc.mmRetina = [arcLength_p1p2(0,deg2rad(axes.opticDisc.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.opticDisc.degRetina(2))), 0];



end

