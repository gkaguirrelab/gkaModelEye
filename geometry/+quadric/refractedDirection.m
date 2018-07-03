function Rr = refractedDirection(R,N,nRel)
%
% Inputs:
%   R                     - 3x2 matrix that specifies the incident ray as a
%                           unit vector of the form [p; d]:
%                               R = p + t*u,
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%   N                     - 3x2 matrix that specifies the surface normal as
%                           a unit vector of the form [p; d]:
%                               N = p + t*u,
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%   nRel                  - The relative index of refraction n / n', where
%                           n is the index of refraction of the medium of
%                           the incident ray, and n' is the index of
%                           refraction of the medium of the surface that
%                           the ray has intersected.
%
% Outputs:
%   Rr              - 3x2 matrix that specifies the refracted ray as 
%                           a unit vector of the form [p; d]:
%                               Rr = p + t*u,
%                           where p is vector origin, d is the direction
%                           expressed as a unit steo, and t has an
%                           obligatory value of unity.
%
% Examples:
%{
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    S = quadric.scale(S,[10 10 10]);
    S = quadric.translate(S,[0; 0; 22]);
    p = [0;0;0];
    u = [0;tand(17.309724);1];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    [X1,X2] = quadric.intersectRay(S,R);
    N = quadric.surfaceNormal(S,X);
    Rr = quadric.refractedDirection(R,Rr,1.2)
%}

% Pre-allocate the output variables
Rrefract = nan(3,2);

% Obtain the direction vector of the incident ray
Ru = R(:,2);

% Obtain the direction vector of the surface normal
Nu = N(:,2);

% Place the surface intersection point as the origin of the refracted ray
Rr(:,1) = N(:,1);

% Calculate the direction vector of the refracted ray
Rr(:,2) = nRel*Ru + nRel*(-Nu*Ru-sqrt(1+(nRel^2)*((Nu*Ru)^2-1)))*Q;

end

