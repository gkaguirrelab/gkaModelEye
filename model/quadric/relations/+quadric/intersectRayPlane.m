function X = intersectRayPlane(n,p,R)
% Find the coordinates of intersection of a ray with a plane
%
% Syntax:
%  X = quadric.intersectRay(S,R,side,boundingBox,bbTol)
%
% Description:
%   Returns the coordinates of the point of intersection of a ray with a
%   plane, if such an intersection exists. Adapted from code written by
%   Nassim Khaled of Wayne State University.
%
% Inputs:
%   n                     - 3x1 vector that is the normal of the plane.
%   p                     - 3x1 vector that is any point in the plane.
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Outputs:
%   X                     - 3x1 vector that gives the coordinates of the
%                           intersection of the ray with the plane
%

% Intiate X for a premature exit
X = [nan; nan; nan];

u = R(:,2);
w = R(:,1) - p;
D = dot(n,u);
N = -dot(n,w);

if abs(D) < 10^-7        % The segment is parallel to plane
    return
end

%compute the intersection parameter
sX = N / D;
X = R(:,1)+ sX.*u;


end

