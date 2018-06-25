function S = canonicalToMatrix(c)


center = c(1:3);
radii = c(4:6);
angles = deg2rad(c(7:9));
scale = c(10);

D = diag(radii);
sgns = sign(D);
Ainv = eul2rotm(angles)*D;
A = inv(Ainv);
Q = A'*A*scale;
Q = Q.*sgns;

% Restore the sign of the radii


Q(4,4) = -scale;

% Translation
T = eye( 4 );
T( 4, 1:3 ) = center;
sgns = ones(4,4);
sgns(4,1:3) = sign(center);
sgns(1:3,4) = sign(center);

S = T * Q * T' .* sgns;

end