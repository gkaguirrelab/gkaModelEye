function [d, theta] = arcLength(S,p1,p2,surfaceTolerance)
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
%   d                     - Scalars. Arc length distance in word coordinate
%                           units
%   theta                 - Scalar in radians. Angle between p1 and p2
%                           w.r.t. the center of the ellipsoid.
%
% Examples:
%{
%}


% Keep the compiler happy by excluding prohibited calls
coder.extrinsic('warning')

% Handle incomplete input arguments
if nargin==3
    surfaceTolerance=1e-10;
end

% Pre-allocate the output variables
d = nan; theta = nan;

% Could have test here that S corresponds to an ellipsoid
%% ADD THIS

% Optionally test that p1 and p2 are on the surface of the quadric
if ~isempty(surfaceTolerance)
    funcS = quadric.vecToFunc(quadric.matrixToVec(S));
    if abs(funcS(p1(1),p1(2),p1(3))) > surfaceTolerance || abs(funcS(p2(1),p2(2),p2(3))) > surfaceTolerance
        warning('One or more passed points are not on the quadric surface within tolerance (%f)',surfaceTolerance);
        return
    end
end

% Obtain the coordinates of the center of the quadric surface
centerS = quadric.center(S);

% The points p1, p2, and c1 define a plane. Obtain the radii of the ellipse
% that is formed by the intersection of S and the plane
[semiMajor,semiMinor,centerX,centerY,centerZ] = planeEllipsoidIntersect(S,p1,p2,centerS);

% Obtain the angles theta1 and theta2. Theta1 is the angle between the apex of the
% ellipsoid and p1, w.r.t. the center of the ellipsoid.

% normalized vectors
p1 = quadric.normalizeVector(p1 - centerS);
p2 = quadric.normalizeVector(p2 - centerS);

% compute angle
theta = acos(dot(p1, p2, 2));

% setup the elliptic integral
ellipticIntegral=@(theta) sqrt(1-sqrt(1-eye.posteriorChamber.radii(2).^2/eye.posteriorChamber.radii(1).^2)^2.*(sin(theta)).^2);
ellipticIntegral_p1p3=@(theta) sqrt(1-sqrt(1-eye.posteriorChamber.radii(3).^2/eye.posteriorChamber.radii(1).^2)^2.*(sin(theta)).^2);
arcLength_p1p2 = @(theta1,theta2) eye.posteriorChamber.radii(1).*integral(ellipticIntegral_p1p2, theta1, theta2);
arcLength_p1p3 = @(theta1,theta2) eye.posteriorChamber.radii(1).*integral(ellipticIntegral_p1p3, theta1, theta2);

% For the calculation, the first theta value is zero, as we are
% calculating distance from the posterior chamber apex (i.e., the
% intersection of the optical axis with the retina).
axes.visual.mmRetina = [arcLength_p1p2(0,deg2rad(axes.visual.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.visual.degRetina(2))), 0];
axes.opticDisc.mmRetina = [arcLength_p1p2(0,deg2rad(axes.opticDisc.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.opticDisc.degRetina(2))), 0];



end

