function [S] = varsToMatrix(A, B, C, D, E, F, G, H, I, K)
% Assembles a quadric matrix from a set of scalar variables
%
% Syntax:
%  [S] = quadric.varsToMatrix(A, B, C, D, E, F, G, H, I, K)
%
% Description:
%   Convert between two forms of expression of the quadric surface.
%
%   The implicit form of a second-order (quadric) surface is:
%       S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%                   2Dxy + 2Exz + 2Fyz +
%                   2Gx + 2Hy + 2Iz + K == 0
%
% 	Note that the order of the cross-terms is xy, xz, yz
%
%   Matrix form:
%       [A D E G;
%        D B F H;
%        E F C I;
%        G H I K]
%
%   Vector form:
%       [A B C D E F G H I K]
%
% Inputs:
%   A, B, C, ...          - The set of quadric parameters as scalars.
%
% Outputs:
%   S                     - 4x4 matrix of the quadric surface.
%

S = [A D E G; D B F H; E F C I; G H I K];

end