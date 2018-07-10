function [A, B, C, D, E, F, G, H, I, K] = matrixToVars(S)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%
%  S = [A D E G; D B F H; E F C I; G H I K];

A = S(1,1);
B = S(2,2);
C = S(3,3);
D = S(2,1);
E = S(3,1);
F = S(3,2);
G = S(4,1);
H = S(4,2);
I = S(4,3);
K = S(4,4);

end