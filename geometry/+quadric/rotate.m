function S = rotate( S, r )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Store the initial scaler value for the quadric
k = S(end,end);

% Find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% Form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% Translate to the center
S = T * S * T';

% Construct the rotation matrix from the Euler angles
R = eul2rotm(deg2rad(r));

% Apply the rotation
S(1:3,1:3) = R*S(1:3,1:3)*R';

% Translate back to original position
S = T * S * T';

% Return to original scale
S = quadric.normalize(S) .* abs(k);

end

