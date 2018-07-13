function X = center( S )
% Center of a central quadric (quadratic surface). When reducing the
% surface to its canonical form, a translation of the coordinate system as
% given by the coordinates of the center must be made.
%
% 0 = f(x,y,z) = a*x^2 + b*y^2 + c*z^2
%              + 2*f*y*z + 2*g*x*z + 2*h*x*y
%              + 2*p*x + 2*q*y + 2*r*z + d
%
% A center of a quadric surface is a point P with the property that any
% line through P
% * determines a chord of the surface whose midpoint is P, or
% * has no point in common with the surface, or
% * lies entirely in the surface.
%
% Input arguments:
% u:
%    parameters of quadric as [x.^2, y.^2, z.^2, x.*y, x.*z, y.*z, x, y, z, 1]
%
% Output arguments:
% x, y, z:
%    coordinates of center

% Copyright 2012 Levente Hunyadi

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

[A, B, C, D, E, F, G, H, I, ~] = quadric.matrixToVars(S);

% E = ...
%   [ a, h, g, p ...
%   ; h, b, f, q ...
%   ; g, f, c, r ...
%   ; p, q, r, d ...
%   ];

P = ...  % minor of E belonging to p
    [ D, E, G ...
    ; B, F, H ...
    ; F, C, I ...
    ];

Q = ...  % minor of E belonging to q
    [ A, E, G ...
    ; D, F, H ...
    ; E, C, I ...
    ];

R = ...  % minor of E belonging to r
    [ A, D, G ...
    ; D, B, H ...
    ; E, F, I ...
    ];

W = ...  % minor of E belonging to d
    [ A, D, E ...
    ; D, B, F ...
    ; E, F, C ...
    ];

d = det(W);
x = -det(P) / d;
y =  det(Q) / d;
z = -det(R) / d;
X = [x; y; z];



end

