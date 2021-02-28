function geodetic = cartToEllipsoidalGeo( X, S )
% Converts Cartesian to ellipsoidal, Jacobian geodetic coordinates
%
% Syntax:
%  geodetic = quadric.cartToEllipsoidalGeo( X, S )
%
% Description:
%   Converts from Cartesian (x, y, z) coordinates on the ellipsoidal
%   surface to ellipsoidal, orthogonal geodetic (beta, gamma, elevation)
%   coordinates.
%
%   The routine takes a coordinate (X) and a quadric (S). A geodetic
%   coordinate is returned of the form beta (latitude), omega (longitude),
%   and elevation=0. The geodetic coordinates are with reference to a
%   centered, non-rotated ellipsoid, with the axes arranged in a standard
%   form such that they are in descending order of length (i.e., semi-axes
%   ordered a => b => c). The variable S can be supplied in either vector
%   or matrix form. The quadric (and associated point X) is translated to
%   place the center at the origin, and rotated to be axis-aligned. The
%   semi-axes of the quadric and point are re-ordered to match the standard
%   form and the geodetic coordinates computed.
%
%   The distinction between parametric and orthogonal geodetic coordinates
%   on the ellipsoidal surface is discussed here:
%
%       https://geographiclib.sourceforge.io/html/triaxial.html
%
%   The calculation is based upon a result presented in:
%
%       Bektas, Sebahattin. "Geodetic computations on triaxial ellipsoid."
%       International Journal of Mining Science (IJMS) 1.1 (2015): 25-34.
%
% Inputs:
%   X                     - 3x1 vector containing the [x, y, z] coordinates
%                           of the point.
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   geodetic              - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%
% Examples:
%{
    %% Confirm the invertibility of the transform
    % Define an ellipsoidal surface
    S = quadric.scale(quadric.unitSphere,[4,2,5]);
    % Pick a point on the surface
    G = [90; -65; 0];
    X = quadric.ellipsoidalGeoToCart( G , S );
    Gprime = quadric.cartToEllipsoidalGeo( X, S );
    assert(max(abs(G-Gprime)) < 1e-3);
%}
%{
    %% Confirm the invertibility of the transform across quadrants
    % Define an ellipsoidal surface
    S = quadric.scale(quadric.unitSphere,[4.7,3,16]);
    % Find a point on the surface by intersecting a ray
    p = [0;0;0];
    u = [1;tand(15);tand(15)];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRayQuadric(S,R);
    for p1=[-1 1]
        for p2=[-1 1]
            for p3=[-1 1]
                quadrant = [p1; p2; p3];
                geodetic = quadric.cartToEllipsoidalGeo( X.*quadrant, S );
                Xprime = quadric.ellipsoidalGeoToCart( geodetic, S);
                assert(max(abs((X.*quadrant)-Xprime)) < 1e-6);
            end
        end
    end
%}


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Translate the quadric so that it is centered at the origin. Translate
% the point accordingly.
origCenter = quadric.center(S);
S = quadric.translate(S,-origCenter);
X = X-origCenter;

% Rotate the quadric so that it is aligned with the cardinal axes and in
% the cardinal orientation.
[S, rotMat] = quadric.alignAxes(S);

% Subject the Cartesian point to the corresponding rotation.
X = (X'*rotMat)';

% The geodetic is undefined exactly at the umbilical points on the
% ellipsoidal surface. Here, we detect the presence of zeros in the
% cartesian coordinate and shift these to be slightly non-zero.
X(X(1:3)==0) = 1e-6;

% Obtain the radii of the quadric surface and distribute the values. We
% adopt the canonical order of a => b => c, but the order returned after
% alignment of the axes is a <= b <= c. This is why c is mapped to the
% first value in the radii, and why the Cartesian coordinate is assembled
% as [z y x]
radii = quadric.radii(S);
a=radii(3);b=radii(2);c=radii(1);
x=X(3);y=X(2);z=X(1);

% Save a store the sign of the y coordinate. We perform the computation
% with abs(y) and thus for omega values in the range 0-180. We then apply
% the sign of y to the omega value after the computation.
ySign = sign(y);
y = abs(y);

% Examine the signs of the Cartesian coordinate to supply the first guess
betaLast = 45.*sign(z);
omegaLast = 45 + 90.*(sign(x)==-1);

% Constants for the iterative solution
nIterations= 20; % number of loops to refine the estimate
ro=180/pi; % convert degrees to radians

% Define a quantity used repeatedly for the calculation
p = sqrt((a^2-b^2)/(a^2-c^2));

% Refine the estimates of beta and omega
for ii=1:nIterations
    beta = rad2deg(atan2( b*z*sin(omegaLast/ro),c*y*sqrt(1-(p^2*cos(omegaLast/ro)^2)) ));
    omega = rad2deg(atan2( a*y*sqrt((p^2-1)*sin(betaLast/ro)^2+1),b*x*cos(betaLast/ro) ));
    betaLast = beta;
    omegaLast = omega;
end

% Mandatory elevation of zero for a point on the surface
elevation = 0;

% Apply the stored sign of the y coordinate to omega
omega = omega.*ySign;

% Assemble the geodetic to return
geodetic=[beta; omega; elevation];

end