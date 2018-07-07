function S = alignAxes(S)

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Store the original scale
Sscale = S(4,4);

% find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% translate to the center
Q = T * S * transpose(T);

% Store the translated scale
Qscale = Q(4,4);

% solve the eigenproblem
[evecs,~] = svd(Q(1:3,1:3) / -Qscale);
Q(1:3,1:3) = evecs'*Q(1:3,1:3)*evecs;

% Restore the Q scale
Q = (Q./Q(4,4)).*Qscale;

% Restore the translation
T = eye( 4 );
T( 4, 1:3 ) = -center';
S = T * Q * transpose(T);

% Restore the original scale
S = (S./S(4,4)).*Sscale;

end