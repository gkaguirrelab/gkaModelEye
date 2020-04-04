function R = normalizeRay( R )
% Normalizes a ray to be a unit vector
%
% Syntax:
%  R = quadric.normalizeRay( R )
%
% Description:
%   Just what it says on the tin.
%
% Inputs:
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Outputs:
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%
% Examples:
%{
    p=[0;0;0];
    u=[4;1;1];
    R=[p,u];
    R=quadric.normalizeRay(R);
%}

% Obtain the direction vector of ray
d = R(:,2);

% Normalize the direction vector
u = d./sqrt(sum(d.^2));

% Assemble the normalized ray
R(:,2) = u;
    
end

