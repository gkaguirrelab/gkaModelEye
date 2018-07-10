function r = radii(S)
% The transparent form of a quadric are the parameters:
%   sx, sy, sz, alpha, beta, gamma, cx, cy, cz, scale
%
% where [sx, sy, sz] are the semi-axes, [alpha, beta, gamma] are the angles
% in degrees, and [cx, cy, cz] is the center.


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% translate to the center
Q = T * S * transpose(T);

% solve the eigenproblem
[~,evals] = svd(Q( 1:3, 1:3 ) / -Q( 4, 4 ));
r = 1./diag(sqrt(evals));
sgns = sign( diag( evals ) );
r = r .* sgns;

r=flipud(r);

end

