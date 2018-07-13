function r = radii(S)


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

% The radius values are at this stage in canonical order (smallest to
% largest). Now re-order the valus so that the correspond to the actual
% x,y,z dimensions of the quadric.
axisOrder = quadric.axisOrder(S);
r = r(axisOrder);

end

