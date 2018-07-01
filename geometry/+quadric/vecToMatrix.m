function S = vecToMatrix(v)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%

A = v(1);
B = v(2);
C = v(3);
D = v(4);
E = v(5);
F = v(6);
G = v(7);
H = v(8);
I = v(9);
K = v(10);

S = [A D E G; D B F H; E F C I; G H I K];

end