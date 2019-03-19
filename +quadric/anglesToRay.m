function R = anglesToRay(p, angle_xy, angle_xz )
% Converts angles with respect to the x axis into a ray
%
% Syntax:
%  R = quadric.anglesToRay( p, angle_xy, angle_xz )
%
% Description:
%   Just what it says on the tin.
%
% Inputs:
%   p                     - 3x1 vector that specifies the ray origin
%   angle_xy, angle_xz    - Scalars. Angles in degrees between the rays
%                           projected on the xy and xz planes.
%
% Outputs:
%   R                     - 3x2 matrix that specifies the refracted ray as 
%                           a unit vector of the form [p; u]:
%                               R = p + t*u,
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%
% Examples:
%{
    % Confirm invertibility
    for angle_xy = [5,-5,-95,95];
        for angle_xz = [5,-5,-95,95];
            p = [-3.7;1;-0.2];
            R = quadric.anglesToRay( p, angle_xy, angle_xz );
            [Rangle_xy, Rangle_xz, Rp] = quadric.rayToAngles( R );
            assert(abs(Rangle_xy-angle_xy)<1e6)
            assert(abs(Rangle_xz-angle_xz)<1e6)
        end
    end
%}

m = 1;
if abs(angle_xy)>90 || abs(angle_xz)>90
    m = -1;
end

d = [m; m*tan(deg2rad(angle_xy)); m*tan(deg2rad(angle_xz))];
R = quadric.normalizeRay([p, d]);

end

