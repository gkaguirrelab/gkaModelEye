function X = center( S )
% Returns the coordinates of the center of the quadric surface
%
% Syntax:
%  X = quadric.center( S )
%
% Description:
%   Given a quadric surface (S) the routine returns the [x; y; z]
%   coordinates of the center of the surface.
%
%   Adapted from a routine written by Levente Hunyadi, who noted that the
%   center of a quadric surface has the property that any line through it:
%     - determines a chord of the surface whose midpoint is P, or
%     - has no point in common with the surface, or
%     - lies entirely in the surface.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   X                     - 3x1 vector containing the [x; y; z] coordinates
%                           of the point.
%
% Examples:
%{
    S = quadric.scale(quadric.unitSphere,[5 3 4]);
    quadric.center(S);
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Decompose the matrix into individual variables
[A, B, C, D, E, F, G, H, I, ~] = quadric.matrixToVars(S);

% Code below taken directly from Levente Hunyadi

% E = ...
%   [ A D E G ...
%   ; D B F H ...
%   ; E F C I ...
%   ; G H I K ...
%   ];

P = ...  % minor of E belonging to g
    [ D, E, G ...
    ; B, F, H ...
    ; F, C, I ...
    ];

Q = ...  % minor of E belonging to h
    [ A, E, G ...
    ; D, F, H ...
    ; E, C, I ...
    ];

R = ...  % minor of E belonging to i
    [ A, D, G ...
    ; D, B, H ...
    ; E, F, I ...
    ];

W = ...  % minor of E belonging to k
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

