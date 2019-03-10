function [angle_xy, angle_xz, p] = rayToAngles( R )
% Converts angles with respect to the x axis into a ray
%
% Syntax:
%  [angle_xy, angle_xz, p] = quadric.rayToAngles( R )
%
% Description:
%   Just what it says on the tin.
%
% Inputs:
%   R                     - 3x2 matrix that specify a vector of the form
%                           [p; d]:
%                               R = p + d,
%                           where p is vector origin, d is the direction.
%
% Outputs:
%   angle_xy, angle_xz    - Scalars. Angles in degrees between the rays
%                           projected on the xy and xz planes.
%   p                     - 3x1 vector that specifies the ray origin
%
% Examples:
%{
    p = [0;0;0];
    u = [1;0;0];
    R1 = quadric.normalizeRay([p, u]);
    p = [0;0;0];
    u = [1;tand(15);tand(-7)];
    R2 = quadric.normalizeRay([p, u]);
    [angle, angle_xy, angle_xz] = quadric.angleRays( R1, R2 )
%}
R = quadric.normalizeRay(R);
p = R(:,1);
u = R(:,2);
angle_xy = rad2deg(atan2(u(2),u(1)));
angle_xz = rad2deg(atan2(u(3),u(1)));
end

