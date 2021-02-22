function distance = distancePointRay(p,R)
% Find the minimum distance between a point and a ray
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
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Outputs:
%   distance              - Scalar distance between the two rays at their
%                           point of closest approach
%


t = dot(R(:,2),(p-R(:,1)))/dot(R(:,2),R(:,2));
d = p-(R(:,1)+t.*R(:,2));
      
% Obtain the Euclidean distance in the 3 dimensions.
distance = norm(d);

end