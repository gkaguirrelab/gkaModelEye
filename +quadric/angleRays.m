function [angle, angle_p1p2, angle_p1p3] = angleRays( R1, R2 )
% Returns the angle between two rays
%
% Syntax:
%  R = quadric.normalizeRay( R )
%
% Description:
%   Just what it says on the tin.
%
% Inputs:
%   R1, R2                - 3x2 matricies that specify a vector of the form
%                           [p; d]:
%                               R = p + d,
%                           where p is vector origin, d is the direction.
%
% Outputs:
%   angle                 - Scalar. Angle in degrees between the rays.
%
% Examples:
%{
    p=[0;0;0];
    u=[4;1;1];
    R=[p,u];
    R=quadric.normalizeRay(R);
%}

angle = rad2deg(atan2(norm(cross(R1(:,2),R2(:,2))), dot(R1(:,2),R2(:,2))));

end

