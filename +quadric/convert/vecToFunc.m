function F = vecToFunc(v)
% Returns a handle to the polymonial function defined by the quadric
%
% Syntax:
%  F = quadric.vecToFunc(v)
%
% Description:
%   A quadric surface is defined by a second order polynomial in three
%   dimensions. This routine returns a handle to the implicit function that
%   defines the quadric surface. The function returns zero when evaluated
%   with an [x, y, z] coordinate on the quadric surface.
%
%	The implicit form of a second-order (quadric) surface:
%       S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%                   2Dxy + 2Exz + 2Fyz +
%                   2Gx + 2Hy + 2Iz + K == 0
%
%   Note that the order of the cross-terms is xy, xz, yz
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%

% If the quadric surface was passed in matrix form, convert to vec
if isequal(size(v),[4 4])
    v = quadric.matrixToVec(v);
end

% Assemble the function
F = @(xx,yy,zz) v(1)*xx.^2 + v(2)*yy.^2 + v(3)*zz.^2 + ...
                2*v(4)*xx.*yy + 2*v(5)*xx.*zz + 2*v(6)*yy.*zz + ...
                2*v(7)*xx + 2*v(8)*yy + 2*v(9)*zz + v(10);

end