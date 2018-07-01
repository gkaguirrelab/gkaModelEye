function [X1, X2] = intersectRay(S,R)
%
% Inputs:
%   S                     - 4x4 quadratic surface matrix
%   R                     - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u,
%                           where p is vector origin, d is the direction
%                           expressed as a unit steo, and t has an
%                           obligatory value of unity.
%
% Outputs:
%   X1, X2                - 3x1 vectors that give the coordinates of the
%                           intersection of the ray with the quadric
%                           surface.
%
% Examples:
%{
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    p = [0; 0; 0];
    u = [0; 1; 1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    [X1, X2] = quadric.intersectRay(S,R);
%}

% Pre-allocate the output variables
X1 = nan(3,1);
X2 = nan(3,1);

% Decompose the ray into homogeneous components
p = [R(:,1); 1];
u = [R(:,2); 1];

% Ray-quadric intersection equation components
%   The variable t of the ray expression p+t*u is given by:
%
%       t = (-b +- sqrt(b^2 - 4ac))/2a
%
% where 
a = u'*S*u;
b = u'*S*p;
c = p'*S*p;

% If the discriminant is less than or equal to zero, then the ray either
% misses or intersects the surface tangentially. We return nans.
if (b^2 - 4*a*c)<=0
    return
end

% If a is close to zero, then obtain the solution using the linear
% approximation:
%   t ~= -c/b
if abs(a)<1e-10
    t(:,1) = c/b;
    t(:,2) = -c/b;
else
% Obtain the quadratic roots.
    t(:,1) = (-b - sqrt(b.^2 - 4*a*c))/(2*a);
    t(:,2) = (-b + sqrt(b.^2 - 4*a*c))/(2*a);
end

% Calculate the coordinates of intersection given t
X1 = p(1:3)+u(1:3)*t(:,1);
X2 = p(1:3)+u(1:3)*t(:,2);

end

