function X = geodeticToCart( geodetic, S )
% Converts geodetic to Cartesian coordinates on an ellipsoidal surface
%
% Syntax:
%  X = quadric.geodeticToCart( geodetic, radii )
%
% Description:
%   Converts from geodetic coordinates (latitude - phi, longitude - lambda,
%   elevation) on the ellipsoidal surface to Cartesian (x, y, z)
%   coordinates.
%
%   The routine takes a coordinate (X) and a quadric (S). A geodetic
%   coordinate is returned of the form latitude (phi), longitude (lambda),
%   and height (distance from the quadric surface). The geodetic
%   coordinates are with reference to a centered, non-rotated ellipsoid,
%   with the axes arranged in a standard form such that they are in
%   descending order of length (i.e., semi-axes ordered a => b => c). The
%   variable S can be supplied in either vector or matrix form. The quadric
%   (and associated point X) is translated to place the center at the
%   origin, and rotated to be axis-aligned. The semi-axes of the quadric
%   and point are re-ordered to match the standard form and the geodetic
%   coordinates computed.
%
%   The operations are modified from a function written by Sebahattin
%   Bektas, (sbektas@omu.edu.tr):
%
%       Bektas, Sebahattin. "Geodetic computations on triaxial ellipsoid."
%       International Journal of Mining Science (IJMS) 1.1 (2015): 25-34.
%
% Inputs:
%   geodetic              - 3x1 vector that provides the geodetic
%                           coordinates latitude, longitude, and elevation
%                           in units of degrees and Cartesian distance. The
%                           latitude is defined over the range -90:90, and
%                           the longitude over the range -180:180.
%                           Elevation takes a value of zero for a point
%                           that is on the surface of ellipsoid.
%   S                     - 1x10 vector or 4x4 matrix specifyin the quadric
%                           surface.
%
% Outputs:
%   X                     - 3x1 vector containing the [x, y, z] coordinates
%                           of the point.
%
% Examples:
%{
    %% Confirm the invertibility of the transform
    % Define an ellipsoidal surface
    S = quadric.scale(quadric.unitSphere,[2,4,5]);
    % Find a point on the surface by intersecting a ray
    p = [0;0;0];
    u = [1;tand(15);tand(-15)];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
    geodetic = quadric.cartToGeodetic( X, quadric.radii(S) );
    Xprime = quadric.geodeticToCart( geodetic, quadric.radii(S) );
    assert(max(abs(X-Xprime)) < 1e-6);
%}


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Translate the quadric so that it is centered at the origin. Store the
% translation values to be applied to the Cartesian coordinate later.
origCenter = quadric.center(S);
S = quadric.translate(S,-origCenter);

% Rotate the quadric so that it is aligned with the cardinal axes.
% Store the rotation matrix to be applied to the Cartesian coordinate later
[S, rotMat] = quadric.alignAxes(S);

% Obtain the radii of the quadric surface and distribute the values.
% Bektas' code expected the radii to be in the order a => b => c, but the
% order returned after alignment of the axes is a <= b <= c. This is why c
% is mapped to the first value in the radii, and why the Cartesian
% coordinate is later assembled as [z y x]
radii = quadric.radii(S);
a=radii(3);b=radii(2);c=radii(1);

% This next block contains essentially unedited code from Bektas.
ro=180/pi; % convert degrees to radians
phi=geodetic(1);lambda=geodetic(2);height=geodetic(3);

ex2=(a^2-c^2)/a^2;
ee2=(a^2-b^2)/a^2;
V=a/sqrt(1-ex2*sin(phi/ro)^2-ee2*cos(phi/ro)^2*sin(lambda/ro)^2);

x=(V+height)*cos(lambda/ro)*cos(phi/ro);
y=(V*(1-ee2)+height)*sin(lambda/ro)*cos(phi/ro);
z=(V*(1-ex2)+height)*sin(phi/ro);

X=[z; y; x];

% Now counter-rotate and then counter-translate the X coordinate
X = (X'*rotMat')';
X = X+origCenter;

end
