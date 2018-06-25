function v = matrixToPolynomial(S)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%

v(1) = S(1,1);
v(2) = S(2,2);
v(3) = S(3,3);
v(4) = S(1,2);
v(5) = S(1,3);
v(6) = S(2,3);
v(7) = S(1,4);
v(8) = S(2,4);
v(9) = S(3,4);
v(10) = S(4,4);

end