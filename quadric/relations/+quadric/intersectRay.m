function [X1, X2] = intersectRay(S,R,side,boundingBox,bbTol)
% Find the coordinates of intersection of a ray with a quadric surface
%
% Syntax:
%  [X1, X2] = quadric.intersectRay(S,R,side,boundingBox,bbTol)
%
% Description:
%   Returns the coordinates of the points of intersection of a ray with a
%   quadric surface.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%   side                  - Scalar, taking the value of -1 or 1. Defines
%                           whether the first (-1) or second (+1)
%                           intersection point should be returned as X1.
%                           For an ellipsoid, setting side=-1 causes the
%                           ray to intersect the convex surface, and
%                           setting side=+1 causes the ray to intersect the
%                           concave surface. For a two-sheeted hyperboloid,
%                           with the ray traveling down the open axis,
%                           setting side = -1 means that the first point of
%                           intersection will be within the interior of one
%                           of the hyperboloid sheets, and thus an
%                           intersection with a concave surface.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates of an intersection are reported.
%                           Intersection points outside of these bounds are
%                           returned as nans. If not defined no bounds will
%                           be enforced.
%   bbTol                 - Scalar. Defines the tolerance within which the
%                           intersection must be within the boundingBox.
%                           Default value is 0.1. Handles the situation in
%                           which the intersection is right on the boundary
%                           but is numerically outside.
%
% Outputs:
%   X1, X2                - 3x1 vectors that give the coordinates of the
%                           intersection of the ray with the quadric
%                           surface.
%
% Examples:
%{
    % Unit sphere, diagonal ray starting from the origin
    S = quadric.unitSphere;
    p = [0; 0; 0];
    u = [1; 1; 1];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRay(S,R);
%}
%{
    % Unit sphere, axis-aligned ray starting from the origin
    S = quadric.unitSphere;
    p = [0; 0; 0];
    u = [0; 0; 1];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRay(S,R);
%}
%{
    % Unit sphere, non-intersecting ray
    S = quadric.unitSphere;
    p = [3; 3; 3];
    u = [0; 0; 1];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRay(S,R);
%}
%{
    % Scaled, translated sphere. Ray starts from sphere center.
    S = quadric.unitSphere;
    S = quadric.scale(S,[2 2 2]);
    S = quadric.translate(S,[0; 1; 1]);
    p = [0;1;1];
    u = [0;tand(17.309724);1];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRay(S,R);
%}
%{
    % Centered ellipsoid. Ray starts from center.
    S = quadric.unitSphere;
    S = quadric.scale(S,[5 4 3]);
    p = [0;0;0];
    u = [1;tand(15);0];
    R = quadric.normalizeRay([p, u]);
    X = quadric.intersectRay(S,R);
%}

% Handle incomplete input arguments
if nargin==2
    side=1;
    boundingBox = [];
    bbTol = 1e-2;
end

if nargin==3
    boundingBox = [];
    bbTol = 1e-2;
end

if nargin==4
    bbTol = 1e-2;
end

if nargin==5
    bbTol = 1e-2;
end

% Pre-allocate the output variables
X1 = nan(3,1);
X2 = nan(3,1);

% Decompose the ray into homogeneous components
p = R(:,1);
u = R(:,2);
px=p(1); py=p(2); pz=p(3);
ux=u(1); uy=u(2); uz=u(3);

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Decompose the quadric matrix into individual variables
[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVars(S);

% Code to generate the ray-quadric intersection equations
%{
    syms A B C D E F G H I K
    syms px py pz ux uy uz t
    % Quadric surface matrix
    S = quadric.varsToMatrix(A, B, C, D, E, F, G, H, I, K);
    % Ray matrix
    p = [px; py; pz];
    u = [ux; uy; uz];
    R = p + t * u;
    % Ray-quadric intersection
    eq = transpose([R; 1]) * S * [R; 1] == 0;
    % Solve for t
    tVals = solve(eq,t);
    % Express as a matlab function
    tVals_fxnHandle = matlabFunction(tVals(1));
%}

% Group the tVals solution into three components. The overall solution is
% of the form
%   t = (-alpha +- sqrt(beta))/gamma
alpha = G.*ux+H.*uy+I.*uz+A.*px.*ux+B.*py.*uy+D.*px.*uy+D.*ux.*py+C.*pz.*uz+E.*px.*uz+E.*ux.*pz+F.*py.*uz+F.*uy.*pz;
beta = G.^2.*ux.^2+H.^2.*uy.^2+I.^2.*uz.^2+D.^2.*px.^2.*uy.^2+D.^2.*ux.^2.*py.^2+E.^2.*px.^2.*uz.^2+E.^2.*ux.^2.*pz.^2+F.^2.*py.^2.*uz.^2+F.^2.*uy.^2.*pz.^2-A.*K.*ux.^2-B.*K.*uy.^2-C.*K.*uz.^2-A.*H.*ux.^2.*py.*2.0-B.*G.*px.*uy.^2.*2.0-C.*G.*px.*uz.^2.*2.0-A.*I.*ux.^2.*pz.*2.0+D.*G.*ux.^2.*py.*2.0+D.*H.*px.*uy.^2.*2.0-C.*H.*py.*uz.^2.*2.0-B.*I.*uy.^2.*pz.*2.0+E.*G.*ux.^2.*pz.*2.0+E.*I.*px.*uz.^2.*2.0+F.*H.*uy.^2.*pz.*2.0+F.*I.*py.*uz.^2.*2.0-A.*B.*px.^2.*uy.^2-A.*B.*ux.^2.*py.^2-A.*C.*px.^2.*uz.^2-A.*C.*ux.^2.*pz.^2-B.*C.*py.^2.*uz.^2-B.*C.*uy.^2.*pz.^2+G.*H.*ux.*uy.*2.0-D.*K.*ux.*uy.*2.0-E.*K.*ux.*uz.*2.0+G.*I.*ux.*uz.*2.0-F.*K.*uy.*uz.*2.0+H.*I.*uy.*uz.*2.0+A.*H.*px.*ux.*uy.*2.0+B.*G.*ux.*py.*uy.*2.0+A.*I.*px.*ux.*uz.*2.0-D.*G.*px.*ux.*uy.*2.0+C.*G.*ux.*pz.*uz.*2.0-D.*H.*ux.*py.*uy.*2.0-E.*G.*px.*ux.*uz.*2.0+B.*I.*py.*uy.*uz.*2.0+C.*H.*uy.*pz.*uz.*2.0-F.*G.*px.*uy.*uz.*4.0+F.*G.*ux.*py.*uz.*2.0+F.*G.*ux.*uy.*pz.*2.0+E.*H.*px.*uy.*uz.*2.0-E.*H.*ux.*py.*uz.*4.0+E.*H.*ux.*uy.*pz.*2.0+D.*I.*px.*uy.*uz.*2.0+D.*I.*ux.*py.*uz.*2.0-D.*I.*ux.*uy.*pz.*4.0-F.*H.*py.*uy.*uz.*2.0-E.*I.*ux.*pz.*uz.*2.0-F.*I.*uy.*pz.*uz.*2.0-A.*F.*ux.^2.*py.*pz.*2.0-A.*F.*px.^2.*uy.*uz.*2.0-B.*E.*px.*uy.^2.*pz.*2.0-C.*D.*px.*py.*uz.^2.*2.0-B.*E.*ux.*py.^2.*uz.*2.0-C.*D.*ux.*uy.*pz.^2.*2.0+D.*E.*ux.^2.*py.*pz.*2.0+D.*E.*px.^2.*uy.*uz.*2.0+D.*F.*px.*uy.^2.*pz.*2.0+D.*F.*ux.*py.^2.*uz.*2.0+E.*F.*px.*py.*uz.^2.*2.0+E.*F.*ux.*uy.*pz.^2.*2.0-D.^2.*px.*ux.*py.*uy.*2.0-E.^2.*px.*ux.*pz.*uz.*2.0-F.^2.*py.*uy.*pz.*uz.*2.0+A.*B.*px.*ux.*py.*uy.*2.0+A.*C.*px.*ux.*pz.*uz.*2.0+A.*F.*px.*ux.*py.*uz.*2.0+A.*F.*px.*ux.*uy.*pz.*2.0+B.*C.*py.*uy.*pz.*uz.*2.0+B.*E.*px.*py.*uy.*uz.*2.0+B.*E.*ux.*py.*uy.*pz.*2.0+C.*D.*px.*uy.*pz.*uz.*2.0+C.*D.*ux.*py.*pz.*uz.*2.0-D.*E.*px.*ux.*py.*uz.*2.0-D.*E.*px.*ux.*uy.*pz.*2.0-D.*F.*px.*py.*uy.*uz.*2.0-D.*F.*ux.*py.*uy.*pz.*2.0-E.*F.*px.*uy.*pz.*uz.*2.0-E.*F.*ux.*py.*pz.*uz.*2.0;
gamma = A.*ux.^2+B.*uy.^2+C.*uz.^2+D.*ux.*uy.*2.0+E.*ux.*uz.*2.0+F.*uy.*uz.*2.0;

% Here, beta is the discriminant. If it is less than zero, then the ray
% misses the surface. In this situation, return nans.
if beta<0
	return
end

% The scale of the ray (t) is given by a quadratic equation. Obtain the two
% roots.
t1 = (-alpha - sqrt(beta))/gamma;
t2 = (-alpha + sqrt(beta))/gamma;

% Calculate the two coordinates of intersection. Order these by the side
% variable
if side==-1
    X1 = p(1:3)+u(1:3)*t1;
    X2 = p(1:3)+u(1:3)*t2;
end
if side==1
    X1 = p(1:3)+u(1:3)*t2;
    X2 = p(1:3)+u(1:3)*t1;
end

% Set to nan coordinates that are not within the bounding box
if ~isempty(boundingBox)
    if ~(X1(1)>=(boundingBox(1)-bbTol) && X1(1)<=(boundingBox(2)+bbTol) && ...
            X1(2)>=(boundingBox(3)-bbTol) && X1(2)<=(boundingBox(4)+bbTol) && ...
            X1(3)>=(boundingBox(5)-bbTol) && X1(3)<=(boundingBox(6)+bbTol))
        X1 = nan(3,1);
    end
    if ~(X2(1)>=(boundingBox(1)-bbTol) && X2(1)<=(boundingBox(2)+bbTol) && ...
            X2(2)>=(boundingBox(3)-bbTol) && X2(2)<=(boundingBox(4)+bbTol) && ...
            X2(3)>=(boundingBox(5)-bbTol) && X2(3)<=(boundingBox(6)+bbTol))
        X2 = nan(3,1);
    end
end

end

