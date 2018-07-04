function [X1, X2] = intersectRay(S,R,side,boundingBox)
%
%
% Syntax:
%  [X1, X2] = intersectRay(S,R)
%
% Description:
%   Returns the coordinates of the points of intersection of a ray with a
%   quadric surface.
%
% Inputs:
%   S                     - 4x4 quadric surface matrix
%   R                     - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%
% Outputs:
%   X1, X2                - 3x1 vectors that give the coordinates of the
%                           intersection of the ray with the quadric
%                           surface.
%
% Examples:
%{
    % Unit sphere, diagonal ray starting from the origin
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    p = [0; 0; 0];
    u = [1; 1; 1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
%}
%{
    % Unit sphere, axis-aligned ray starting from the origin
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    p = [0; 0; 0];
    u = [0; 0; 1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
%}
%{
    % Unit sphere, non-intersecting ray
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    p = [3; 3; 3];
    u = [0; 0; 1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
%}
%{
    % Scaled, translated sphere. Ray starts from sphere center.
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    S = quadric.scale(S,[2 2 2]);
    S = quadric.translate(S,[0; 1; 1]);
    p = [0;1;1];
    u = [0;tand(17.309724);1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    X = quadric.intersectRay(S,R);
%}

if nargin==2
    side=1;
end

if nargin==3
    boundingBox = [];
end

% Pre-allocate the output variables
X1 = nan(3,1);
X2 = nan(3,1);

% Decompose the ray into homogeneous components
p = R(:,1);
u = R(:,2);
px=p(1); py=p(2); pz=p(3);
ux=u(1); uy=u(2); uz=u(3);

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

% Here, beta is the discriminant. If it is less than or equal to zero, then
% the ray either misses or intersects the surface tangentially. In this
% situation, return nans.
if beta<=0
    return
end

% The scale of the ray (t) is given by a quadratic equation. Obtain the two
% roots.
t(:,1) = (-alpha - sqrt(beta))/gamma;
t(:,2) = (-alpha + sqrt(beta))/gamma;

% Calculate the two coordinates of intersection. Order these by the side
% variable
if side==1
    X1 = p(1:3)+u(1:3)*t(:,1);
    X2 = p(1:3)+u(1:3)*t(:,2);
end
if side==2
    X1 = p(1:3)+u(1:3)*t(:,2);
    X2 = p(1:3)+u(1:3)*t(:,1);
end

% If a bounding box has been specified, check to see if only one of the
% coordinates is within the box, and if so, report this one as X1
if ~isempty(boundingBox)
    inX1 = false;
    if X1(1)>=boundingBox(1) && X1(1)<=boundingBox(2) && ...
            X1(2)>=boundingBox(3) && X1(2)<=boundingBox(4) && ...
            X1(3)>=boundingBox(5) && X1(3)<=boundingBox(6)
        inX1 = true;
    else
            X1 = nan(3,1);
    end
    inX2 = false;
    if X2(1)>=boundingBox(1) && X2(1)<=boundingBox(2) && ...
            X2(2)>=boundingBox(3) && X2(2)<=boundingBox(4) && ...
            X2(3)>=boundingBox(5) && X2(3)<=boundingBox(6)
        inX2 = true;
    else
            X2 = nan(3,1);
    end
    if inX2 && ~inX1
        tmp = X1;
        X1 = X2;
        X2 = tmp;
    end
end

end

