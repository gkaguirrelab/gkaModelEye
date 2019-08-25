function distance = distancePointRay(p,R)
% Find the minimum distance between a ray and a point
%
% Syntax:
%  distance = quadric.distancePointRay(p,R)
%
% Description:
%   Yep. That's what it does.
%
% Inputs:
%   p                     - 3x1 vector that specifies a point
%   R                     - 3x2 matrix that specifies a vector of the form
%                           [p; d]:
%                               R = p + d,
%                           where p is vector origin, d is the direction.
%
% Outputs:
%   distance              - Scalar distance between the two rays at their
%                           point of closest approach
%

d = norm(cross(R(:,2),p - R(:,1))) ...
    / norm(R(:,1));
      
% Obtain the Euclidean distance in the 3 dimensions.
distance = sqrt(sum(d.^2));

end