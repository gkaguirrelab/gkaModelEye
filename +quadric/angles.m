function angles = angles(S)
% The transparent form of a quadric are the parameters:
%   sx, sy, sz, alpha, beta, gamma, cx, cy, cz, scale
%
% where [sx, sy, sz] are the semi-axes, [alpha, beta, gamma] are the angles
% in degrees, and [cx, cy, cz] is the center.


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% solve the eigenproblem
[evecs,~] = svd(-S( 1:3, 1:3 ) );

% Derive the angles
alfax=atan2(-evecs(3,2),evecs(3,3));
alfay=acos(evecs(3,3)/cos(alfax));
alfaz=atan2(-evecs(2,1),evecs(1,1));
angles=180/pi.*[alfax,alfay,alfaz];

end

