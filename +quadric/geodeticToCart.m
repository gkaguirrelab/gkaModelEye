function X = geodeticToCart( geodetic, radii )
% Converts geodetic to Cartesian coordinates on an ellipsoidal surface
%
% Syntax:
%  X = quadric.geodeticToCart( geodetic, radii )
%
% Description:
%   Converts from geodetic coordinates (latitude, longitude, elevation) on
%   the ellipsoidal surface to Cartesian (x, y, z) coordinates.
%
%   The coordinates are with reference to a centered, aligned ellipsoid.
%   The radii must be provided in the "canonical" order returned by
%   quadric.radii(). The x,y,z coordinates returned are with respect to
%   this order [a <= b <= c].
%
%   The operations are taken from a function written by Sebahattin Bektas,
%   (sbektas@omu.edu.tr) 19 Mayis University, Samsun
%
% Inputs:
%   geodetic              - 3x1 vector that provides the geodetic
%                           coordinates latitude, longitude, and elevation
%                           in units of degrees and Cartesian distance. The
%                           latitude is defined over the range -90:90, and
%                           the longitude over the range -180:180.
%                           Elevation takes a value of zero for a point
%                           that is on the surface of ellipsoid.
%   radii                 - 3x1 vector. Semi-axis lengths of the ellipsoid,
%                           provided in canonical size order (a <= b <= c).
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
    u = [1;tand(15);0];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
    geodetic = quadric.cartToGeodetic( X, quadric.radii(S) );
    Xprime = quadric.geodeticToCart( geodetic, quadric.radii(S) )
    assert(max(abs(X-Xprime)) < 1e-6);
%}


ro=180/pi; % convert degrees to radians

% Distribute the radii to the variables a, b, c. Bektas' original code
% expected the radii in the opposite order. This is why c is mapped to the
% first value in the radii, and why the returned Cartesian coordinate is
% assembled as [z y x]
c=radii(1);b=radii(2);a=radii(3);

fi=geodetic(1);l=geodetic(2);h=geodetic(3);

ex2=(a^2-c^2)/a^2;
ee2=(a^2-b^2)/a^2;
V=a/sqrt(1-ex2*sin(fi/ro)^2-ee2*cos(fi/ro)^2*sin(l/ro)^2);

x=(V+h)*cos(l/ro)*cos(fi/ro);
y=(V*(1-ee2)+h)*sin(l/ro)*cos(fi/ro);
z=(V*(1-ex2)+h)*sin(fi/ro);

X=[z; y; x];

end
