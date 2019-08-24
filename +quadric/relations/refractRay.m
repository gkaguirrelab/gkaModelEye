function Rr = refractRay(R,N,nRel)
% Returns the refracted ray given incident ray, surface normal, and index 
%
% Syntax:
%  Rr = quadric.refractRay(R,N,nRel)
%
% Description:
%   Snell's Law of optical refraction may be expressed in three-dimensions
%   in matrix form. This routine provides the refracted ray after an
%   incident ray encounters a surface with a given surface normal, and
%   possesing a given relative index of refraction. If the ray is found to
%   experience total internal reflection, the routine returns nans for the
%   refracted ray.
%
%   This approach to optical systems is described in:
%
%       Langenbucher, Achim, et al. "Ray tracing through a schematic eye
%       containing second?order (quadric) surfaces using 4× 4 matrix
%       notation." Ophthalmic and Physiological Optics 26.2 (2006):
%       180-188.
%
% Inputs:
%   R                     - 3x2 matrix that specifies the incident ray as a
%                           unit vector of the form [p; u]:
%                               R = p + t*u,
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%   N                     - 3x2 matrix that specifies the surface normal as
%                           a unit vector of the form [p; u]:
%                               N = p + t*u,
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%   nRel                  - The relative index of refraction n / n', where
%                           n is the index of refraction of the medium of
%                           the incident ray, and n' is the index of
%                           refraction of the medium of the surface that
%                           the ray has intersected.
%
% Outputs:
%   Rr                    - 3x2 matrix that specifies the refracted ray as 
%                           a unit vector of the form [p; u]:
%                               Rr = p + t*u,
%                           where p is vector origin, u is the direction
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
    R = quadric.normalizeRay([p, u]);
    n0 = 1;
    n1 = 1.5;
    nRel = n0/n1;
    u=[0;-2;-1];
    N=quadric.normalizeRay([p,u]);
    Rr = quadric.refractRay(R,N,nRel);
    assert(max(abs(Rr(:,2)-[0.629;0.661;0.409]))<0.001);
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

% Compute the full expression in parts to save on computations
Cqu = cross(q,u);
rootTerm = 1-(nRel^2)*dot(Cqu,Cqu);

% Test if the component to be rooted is negative. If so, this reflects
% total internal reflection, and we exit the function with nans for the
% refracted ray.
if rootTerm<0
    return
end

% Complete the calculation
Rr(:,2) = nRel*cross(q,(-Cqu)) - q*sqrt(rootTerm);

end

