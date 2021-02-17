function [angle, angle_xy, angle_xz] = angleRays( R1, R2 )
% Returns the angle in degrees between two rays
%
% Syntax:
%  [angle, angle_xy, angle_xz] = quadric.angleRays( R1, R2 )
%
% Description:
%   Just what it says on the tin.
%
% Inputs:
%   R1, R2                - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Outputs:
%   angle                 - Scalar. Angle in degrees between the rays.
%   angle_xy, angle_xz    - Scalars. Angles in degrees between the rays
%                           projected on the xy and xz planes.
%
% Examples:
%{
    p = [0;0;0];
    u = [1;0;0];
    R1 = quadric.normalizeRay([p, u]);
    p = [0;0;0];
    u = [1;tand(15);tand(-7)];
    R2 = quadric.normalizeRay([p, u]);
    [angle, angle_xy, angle_xz] = quadric.angleRays( R1, R2 );
%}

% Calculate the overall angle between the two vectors
angle = rad2deg(atan2(norm(cross(R1(:,2),R2(:,2))), dot(R1(:,2),R2(:,2))));

% Now obtain the angles as projected on the xy and xz planes
u1_xy = [R1(1,2);R1(2,2);0];
u2_xy = [R2(1,2);R2(2,2);0];
angle_xy = rad2deg(atan2(norm(cross(u1_xy,u2_xy)), dot(u1_xy,u2_xy)));

% Make the angle signed with respect to the R1 vector
angle_xy = angle_xy*sign(dot([0;0;1],cross(u1_xy,u2_xy)));

u1_xz = [R1(1,2);0;R1(3,2)];
u2_xz = [R2(1,2);0;R2(3,2)];
angle_xz = rad2deg(atan2(norm(cross(u1_xz,u2_xz)), dot(u1_xz,u2_xz)));

% Make the angle signed with respect to the R1 vector
angle_xz = angle_xz*sign(dot([0;1;0],cross(u1_xz,u2_xz)));


end

