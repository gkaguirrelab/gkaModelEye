function Rr = refractRay(R,N,nRel)
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
    % Test 1 from http://www.starkeffects.com/snells-law-vector.shtml
    p=[0;0;0];
    u=[1/sqrt(2);0;1/sqrt(2)];
    R=[p,u];
    n0 = 1;
    n1 = 1.5;
    nRel = n0/n1;
    N=[p,[0;0;-1]];
    Rr = quadric.refractRay(R,N,nRel);
    assert(max(abs(Rr(:,2)-[0.471;0;0.882]))<0.001);
%}
%{
    % Test 2 from http://www.starkeffects.com/snells-law-vector.shtml
    p=[0;0;0];
    u=[4;1;1];
    u = u./sqrt(sum(u.^2));
    R=[p,u];
    n0 = 1;
    n1 = 1.5;
    nRel = n0/n1;
    u=[0;-2;-1];
    u = u./sqrt(sum(u.^2));
    N=[p,u];
    Rr = quadric.refractRay(R,N,nRel);
    assert(max(abs(Rr(:,2)-[0.629;0.661;0.409]))<0.001);
%}
%{
    % Elagha 2017 numerical example
    % The paper provides a numerical example in section C. Test that we get
    % the same value.
    % First surface from Elagha 
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    S = quadric.scale(S,[10 10 10]);
    S = quadric.translate(S,[22; 0; 0]);
    p = [0;0;0];
    u = [1;tand(17.309724);0];
    u = u./sqrt(sum(u.^2));
    R = [p, u];
    side = 1;
    X = quadric.intersectRay(S,R,side);
    N = quadric.surfaceNormal(S,X);
    Rr = quadric.refractRay(R,N,1/1.2);
    assert(abs(atan(Rr(2,2)/Rr(1,2))-0.1655)<0.001);
%}

% Pre-allocate the output variable
Rr = nan(3,2);

% Obtain the direction vector of the incident ray
u = R(:,2);

% Obtain the direction vector of the surface normal
q = N(:,2);

% Place the surface intersection point as the origin of the refracted ray
Rr(:,1) = N(:,1);

% Code to generate the 3D expression of Snell's Law
%{
    syms qx qy qz r0x r0y r0z r1x r1y r1z n0 n1
    Q = [qx;qy;qz];
    R0 = [r0x;r0y;r0z];
    R1 = [r1x;r1y;r1z];
    eq = cross(Q,n0*cross(R0,Q)) == cross(Q,n1*cross(R1,Q));

    eq_R = solve(eq,R1);
%}

% Calculate the direction vector of the refracted ray. This equation is
% taken from: http://www.starkeffects.com/snells-law-vector.shtml
Rr(:,2) = nRel*cross(q,(cross(-q,u))) - q*sqrt(1-(nRel^2)*dot(cross(q,u),cross(q,u)));

end

