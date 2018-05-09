% adjointSym3 - recover the adjoint matrix of a 3x3 symmetric matrix
%          - Pierluigi Taddei (pierluigi.taddei@polimi.it)
%
% Usage:   A = adjoint3(M)
%
% Arguments:
%           M - symmetric 3x3 matrix
%           A - the adjoint matrix
% 08.3.2007 : Created
%
function A = adjointSym3(M)
    A = zeros(3,3);
    a = M(1,1); b = M(1,2); d = M(1,3);
                c = M(2,2); e = M(2,3);
                            f = M(3,3);

    A(1,1) = c*f - e*e;
    A(1,2) = -b*f + e*d;
    A(1,3) = b*e - c*d;

    A(2,1) = A(1,2);
    A(2,2) = a*f - d*d;
    A(2,3) = -a*e + b*d;

    A(3,1) = A(1,3);
    A(3,2) = A(2,3);
    A(3,3) = a*c - b*b;
end