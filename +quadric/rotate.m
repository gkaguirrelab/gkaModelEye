function S = rotate( S, r )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here



% If the quadric surface was passed in vector form, convert to matrix
returnVecFlag = false;
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
    returnVecFlag = true;
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

% Construct the rotation matrix from the Euler angles
R = eul2rotm(deg2rad(r));

% Apply the rotation
Q(1:3,1:3) = R'*Q(1:3,1:3)*R;

% Restore the Q scale
Q = (Q./Q(4,4)).*Qscale;

% Restore the translation
T = eye( 4 );
T( 4, 1:3 ) = -center';
S = T * Q * transpose(T);

% Restore the original scale
S = (S./S(4,4)).*Sscale;

% Return a vector if that was the original input
if returnVecFlag
    S = quadric.matrixToVec(S);
end


end

