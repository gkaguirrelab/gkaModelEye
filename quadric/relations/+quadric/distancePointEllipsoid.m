function [distance,Xprojection] = distancePointEllipsoid( X, S )
% Distance and projection of a point onto an ellipsoidal surface
%
% Syntax:
%  [distance,Xprojection] = quadric.distancePointEllipsoid( X, S )
%
% Description:
%   Identifies the closest point on the surface an ellipsoid from a
%   specified point.
%
%   This code is derived from the function "shortest_distance", written by:
%
%       Sebahattin Bektas, 19 Mayis University, Samsun
%       sbektas@omu.edu.tr
%
% Inputs:
%   X                     - 3x1 vector containing the [x, y, z] coordinates
%                           of the point.
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   distance              - Scalar. Distance of the point from the closest
%                           location on the ellipsoidal surface.
%   Xprojection           - 3x1 vector containing the [x, y, z] coordinates
%                           of the closest point on the ellipsoid
%
% Examples:
%{
    % Numeric example from Bektas, Sebahattin. "Orthogonal distance from an
    % ellipsoid." Boletim de Ciências Geodésicas 20.4 (2014): 970-983.
    S = quadric.scale(quadric.unitSphere,[6378388.0000,6378318.0000,6356911.9461]);
    X = [3909863.9271; 3909778.1230; 3170932.5016];
    [distance,Xprojection] = quadric.distancePointEllipsoid( X, S );
    X0projection = [3909251.5547; 3909165.7506; 3170432.5016];
    distance0 = 1e3;
    % Test the result matches Bektas
    assert(norm(Xprojection-X0projection)<1e6)
    assert(abs(distance-distance0)<1e6)
%}

% Confirm that the passed quadric is an ellipsoid
if ~strcmp(quadric.classify(S),'ellipsoid')
    error('distancePointEllipsoid:invalidQuadric','The solution is only valid for ellipsoidal quadrics')
end

% The Bektas routine performs the calculation for a centered, un-rotated
% ellipsoid, in standard orientation. We place the ellipsoid in this
% standard form, and apply the same transformation to the point.

% Translation
T = quadric.center(S);
S = quadric.translate(S,-T);
X = X-T;

% Rotation - step 1
% This places the quadric with the radii ordered short - medium - long with
% the x,y,z axes
[S, evecs] = quadric.alignAxes(S);
X = evecs'*X;

% Rotation - step 2
% Rotate the quadric around the z axis to place the radii in the order that
% the Bektas routine requires: long-medium-short.
S=quadric.rotate( S, [0;90;0] );
X = X([3 2 1]);

% Obtain the radii and place in the a,b,c, variables
radii = quadric.radii(S);
a = radii(1); b = radii(2); c = radii(3);

% Place the coordinate in x,y,z variables
x=X(1);y=X(2);z=X(3);

% Bektas calculations follow
E=sign(a)/a^2;F=sign(b)/b^2;G=sign(c)/c^2;
xo=x;yo=y;zo=z;

% Silence a warning that can occur with the regression step
warnState = warning();
warning('off','MATLAB:nearlySingularMatrix');

for i=1:20
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

% Restore the warn state
warning(warnState);

% The solution
distance = sqrt((x-xo)^2+(y-yo)^2+(z-zo)^2);
Xprojection = [xo; yo; zo];

% Now transform the Xprojection back to the original coordinate space
Xprojection = Xprojection([3 2 1]);
Xprojection = evecs*Xprojection;
Xprojection = Xprojection+T;

end