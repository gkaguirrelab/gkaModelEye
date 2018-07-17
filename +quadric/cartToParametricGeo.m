function geodetic = cartToParametricGeo( X, S )
% Converts Cartesian to parametric geodetic coordinates on an ellipsoid
%
% Syntax:
%   geodetic = quadric.cartToParametricGeo( X, S )
%
% Description:
%   Converts from Cartesian (x, y, z) coordinates on the ellipsoidal
%   surface to parametric geodetic (latitude, longitude, elevation)
%   coordinates.
%
%   The routine takes a coordinate (X) and a quadric (S). A geodetic
%   coordinate is returned of the form phi (latitude), lambda (longitude),
%   and elevation (distance from the quadric surface). The geodetic
%   coordinates are with reference to a centered, non-rotated ellipsoid,
%   with the axes arranged in a standard form such that they are in
%   descending order of length (i.e., semi-axes ordered a => b => c). The
%   variable S can be supplied in either vector or matrix form. The quadric
%   (and associated point X) is translated to place the center at the
%   origin, and rotated to be axis-aligned. The semi-axes of the quadric
%   and point are re-ordered to match the standard form and the geodetic
%   coordinates computed.
%
%   The distinction between parametric and orthogonal geodetic coordinates
%   on the ellipsoidal surface is discussed here:
%
%       https://geographiclib.sourceforge.io/html/triaxial.html
%
%   The operations are modified from a function written by Sebahattin
%   Bektas, (sbektas@omu.edu.tr):
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
%                           coordinates phi, lambda, and elevation in units
%                           of degrees and Cartesian distance. Phi is
%                           defined over the range -90:90, and lambda over
%                           -180:180. Elevation takes a value of zero for a
%                           point that is on the surface of ellipsoid.
%
% Examples:
%{
    %% Confirm the invertibility of the transform
    % Define an ellipsoidal surface
    S = quadric.scale(quadric.unitSphere,[4,2,5]);
    % Find a point on the surface by intersecting a ray
    p = [0;0;0];
    u = [1;tand(15);tand(-15)];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
    geodetic = quadric.cartToParametricGeo( X, S );
    Xprime = quadric.parametricGeoToCart( geodetic, S );
    assert(max(abs(X-Xprime)) < 1e-6);
%}
%{
    %% Confirm the invertibility of the transform across quadrants
    % Define an ellipsoidal surface
    S = quadric.scale(quadric.unitSphere,[4.7,3,16]);
    % Find a point on the surface by intersecting a ray
    p = [0;0;0];
    u = [1;tand(15);tand(15)];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
    for p1=-1:2:1
        for p2=-1:2:1
            for p3=-1:2:1
                quadrant = [p1; p2; p3];
                geodetic = quadric.cartToParametricGeo( X.*quadrant, S );
                Xprime = quadric.parametricGeoToCart( geodetic, S);
                fprintf('X: %d %d %d; geo: %d %d %f; Xp: %d %d %d \n',sign(p1),sign(p2),sign(p3),sign(geodetic(1)),sign(geodetic(2)),geodetic(3),sign(Xprime(1)),sign(Xprime(2)),sign(Xprime(3)));
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

% Obtain the radii of the quadric surface and distribute the values.
% Bektas' code expected the radii to be in the order a => b => c, but the
% order returned after alignment of the axes is a <= b <= c. This is why c
% is mapped to the first value in the radii, and why the Cartesian
% coordinate is assembled as [z y x]
radii = quadric.radii(S);
a=radii(3);b=radii(2);c=radii(1);
x=X(3);y=X(2);z=X(1);

% This next block contains essentially unedited code from Bektas.
% https://www.mathworks.com/matlabcentral/fileexchange/46248-converter-cartesian-coordinate-to-geodetic-coordinate

% Constants
nIterations= 20; % number of loops to refine the estimate
ro=180/pi; % convert degrees to radians
eps=0.0005; % three sholder


ex2=(a^2-c^2)/a^2; ee2=(a^2-b^2)/a^2;
E=1/a^2;F=1/b^2;G=1/c^2;

xo=a*x/sqrt(x^2+y^2+z^2);
yo=b*y/sqrt(x^2+y^2+z^2);
zo=c*z/sqrt(x^2+y^2+z^2);

for ii=1:nIterations
    j11=F*yo-(yo-y)*E;
    j12=(xo-x)*F-E*xo;
    
    j21=G*zo-(zo-z)*E;
    j23=(xo-x)*G-E*xo;
    
    A=[ j11   j12   0
        j21   0   j23
        2*E*xo    2*F*yo  2*G*zo  ];
    
    sa=(xo-x)*F*yo-(yo-y)*E*xo;
    sb=(xo-x)*G*zo-(zo-z)*E*xo;
    se=E*xo^2+F*yo^2+G*zo^2-1;
    Ab=[ sa  sb  se]';
    bil=-A\Ab;
    xo=xo+bil(1);
    yo=yo+bil(2);
    zo=zo+bil(3);
    
    if max(abs(bil))<eps
        break
    end
end

% Modified the Bektas code to use atan2 and thus make all four quadrants of
% the geodesic coordinates invertible for the two angles
phi=ro*atan2(zo*(1-ee2)/(1-ex2),sqrt((1-ee2)^2*xo^2+yo^2));
lambda=ro*atan2(1/(1-ee2)*yo,xo);

% Calculate the elevation of the point w.r.t. the quadric surface
elevation=sign(z-zo)*sign(zo)*sqrt((x-xo)^2+(y-yo)^2+(z-zo)^2);

% Assemble the geodetic to return
geodetic=[phi; lambda; elevation];

end