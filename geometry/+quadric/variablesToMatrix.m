function [S] = variablesToMatrix(A, B, C, D, E, F, G, H, I, K)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%
%  S = [A D E G; D B F H; E F C I; G H I K];

S = [A D E G; D B F H; E F C I; G H I K];

end