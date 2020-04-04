function Rr = reflectRay(R,N)
% Returns the reflected ray given incident ray and surface normal
%
% Syntax:
%  Rr = quadric.reflectRay(R,N)
%
% Description:
%   The angle of reflection equals the angle of incidence!
%
% Inputs:
%   R                     - 3x2 matrix that specifies the incident ray as 
%                           a vector of the form [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%   N                     - 3x2 matrix that specifies the surface normal as 
%                           a vector of the form [p; u], corresponding to
%                               N = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Outputs:
%   Rr                    - 3x2 matrix that specifies the reflected ray as 
%                           a vector of the form [p; u], corresponding to
%                               Rr = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%
% Examples:
%{
    S = quadric.scale(quadric.unitSphere(),[5,8,7]);
    boundingBox = [-10 10 -10 10 -10 10];
    figure
    quadric.plotSurface(S, boundingBox,'k',0.2);
    hold on
    camlight
    R = quadric.anglesToRay([-10; 0; 0], 20, 0 );
    X = quadric.intersectRay(S,R,-1);
    N = quadric.surfaceNormal(S,X,1);
    Rr = quadric.reflectRay(R,N);
    t = abs((X(1)-R(1,1))./R(1,2));
    quadric.plotRay(R,t);
    quadric.plotRay(Rr,2);
%}


% Pre-allocate the output variable
Rr = nan(3,2);

% Clear the nan cases
if any(isnan(R))    
    return
end
if any(isnan(N))
    Rr = R;
    return
end

% Obtain the direction vector of the incident ray
u = R(:,2);

% Obtain the direction vector of the surface normal
q = N(:,2);

% Place the surface intersection point as the origin of the refracted ray
Rr(:,1) = N(:,1);

% Obtain the direction of the reflected ray
Rr(:,2) = u - 2.*(dot(u,q)).*q;

end

