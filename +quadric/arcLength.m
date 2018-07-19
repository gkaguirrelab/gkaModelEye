function [d, theta] = arcLength(S,X1,X2)
% Find the arc distance between two points on a quadric surface
%
% Syntax:
%  d = quadric.arcLength(S,p1,p2)
%
% Description:
%   Returns the arc length distance on the quadric surface between two
%   points. This currently only works for ellipsoids.
%https://geographiclib.sourceforge.io/html/triaxial.html
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

    S = quadric.scale(quadric.unitSphere,[6378172, 6378103, 6356753]);
    % Pick a point on the surface
    G = [0; 0.5; 0];
    X1 = quadric.ellipsoidalGeoToCart( G , S );
    G = [80; 0.5; 0];
    X2 = quadric.ellipsoidalGeoToCart( G , S );
    quadric.arcLength(S,X1,X2);
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Obtain the geodetics for the cartesian coordinates
G1 = quadric.cartToEllipsoidalGeo( X1, S );
G2 = quadric.cartToEllipsoidalGeo( X2, S );

% Identify the angle alpha between the two points, where alpha is the angle
% the geodesic makes with lines of constant omega
alpha = 0;

% Obtain the radii of the quadric surface and distribute the values. We
% adopt the canonical order of a => b => c, but the order returned after
% alignment of the axes is a <= b <= c. This is why c is mapped to the
% first value in the radii.
radii = quadric.radii(quadric.alignAxes(S));
a=radii(3);b=radii(2);c=radii(1);

% Convert degrees to radians
ro=180/pi; 

% Calculate gamma
gamma = (b^2 - c^2)*cos(G1(1)/ro)^2*sin(alpha/ro)^2-(a^2-b^2)*sin(G1(2)/ro)^2*cos(alpha/ro)^2;

ds_dBeta = @(beta) (sqrt(b.^2*sin(beta./ro).^2+c.^2.*cos(beta./ro).^2).*sqrt((b.^2-c.^2)*cos(beta./ro).^2-gamma))./ ...
    sqrt(a.^2-b.^2.*sin(beta./ro).^2-c.^2*cos(beta./ro).^2);
ds_dOmega = @(omega) (sqrt(a.^2.*sin(omega./ro).^2+b.^2*cos(omega./ro).^2).*sqrt((a.^2-b.^2).*sin(omega./ro).^2+gamma))./ ...
    sqrt(a.^2.*sin(omega./ro).^2+b.^2.*cos(omega./ro).^2-c.^2);

s_beta = real(integral(ds_dBeta,G1(1)/ro,G2(1)/ro));
s_omega = real(integral(ds_dOmega,G1(2)/ro,G2(2)/ro));
s = s_beta+s_omega;


end

