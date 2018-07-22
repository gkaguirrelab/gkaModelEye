function X = ellipsoidalGeoToCart( geodetic, S )
% Converts ellipsoidal, Jacobian geodetic to Cartesian coordinates
%
% Syntax:
%  X = quadric.ellipsoidalGeoToCart( geodetic, S )
%
% Description:
%   Converts from the ellipsoidal geodetic coordinates (beta, omega,
%   elevation) on the ellipsoidal surface to Cartesian (x, y, z)
%   coordinates.
%
%   The current implementation functions only for points on the quadric
%   surface. Thus, elevation must be zero.
%
%   The routine takes a geodetic coordinate of the form [beta, omega, 0],
%   and a quadric (S). A coordinate (X) is returned. The geodetic
%   coordinates are with reference to a centered, non-rotated ellipsoid.
%   The variable S can be supplied in either vector or matrix form. The
%   quadric (and associated point X) is translated to place the center at
%   the origin, and rotated to be axis-aligned.
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
%   geodetic              - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation takes a value of zero for a point
%                           that is on the surface of ellipsoid.
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   X                     - 3x1 vector containing the [x, y, z] coordinates
%                           of the point.
%
% Examples:
%{
    % Show an ellipsoidal surface with lines of constant beta and omega
    S = quadric.scale(quadric.unitSphere,[1 1.5 0.5]);
    F = quadric.vecToFunc(quadric.matrixToVec(S));
    figure
    [xx, yy, zz] = meshgrid( linspace(-2,2,100), linspace(-2,2,100), linspace(-2,2,100));
    vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);
    p = patch(vertices);
    p.FaceColor =[.9 .9 .9];
    p.EdgeColor = 'none';
    alpha(.8);
    daspect([1 1 1])
    view([40 30]); 
    axis tight
    axis equal
    camlight
    lighting gouraud
    hold on
    % Plot lines of varying beta
    for beta = -90:10:90
        coords =[];
        for omega = -180:3:180
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        plot3(coords(:,1),coords(:,2),coords(:,3),'-b');
    end
    fprintf('Lines of constant beta in blue\n');
    % Plot lines of varying omega
    for omega = -180:10:180
        coords =[];
        for beta = -90:3:90
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        plot3(coords(:,1),coords(:,2),coords(:,3),'-g');
    end
    fprintf('Lines of constant omega in green\n');
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

% Obtain the radii of the quadric surface and distribute the values. We
% assume the order a => b => c, but the order returned after alignment of
% the axes is a <= b <= c. This is why c is mapped to the first value in
% the radii, and why the Cartesian coordinate is later assembled as [z y x]
radii = quadric.radii(S);

a=radii(3);b=radii(2);c=radii(1);

ro=180/pi; % convert degrees to radians
beta=geodetic(1);omega=geodetic(2);

x=a*cos(omega/ro)*(sqrt(a^2-b^2*sin(beta/ro)^2-c^2*cos(beta/ro)^2)/sqrt(a^2-c^2));
y=b*cos(beta/ro)*sin(omega/ro);
z=c*sin(beta/ro)*(sqrt(a^2*sin(omega/ro)^2+b^2*cos(omega/ro)^2-c^2)/sqrt(a^2-c^2));

% Assemble the coordinate
X=[z; y; x];

% Now counter-rotate and then counter-translate the X coordinate
X = (X'*rotMat')';
X = X+origCenter;

end
