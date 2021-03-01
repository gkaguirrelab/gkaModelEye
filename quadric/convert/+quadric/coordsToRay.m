function R = coordsToRay(pp)
% Converts a pair of coordinates into the ray that connects them
%
% Syntax:
%  R = quadric.coordsToRay(pp)
%
% Description:
%   Given a pair of coordinates, this routine returns the normalized ray
%   that is directed from the first coordinate to the second.
%
% Inputs:
%   pp                    - 3x2 matrix containing two coordinates
%
% Outputs:
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Examples:
%{
%}

R = quadric.normalizeRay([pp(:,1), pp(:,2)-pp(:,1)]);

end

