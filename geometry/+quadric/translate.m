function S = translate( S, t )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Store the initial scaler value for the quadric
k = S(end,end);

% Form the translation matrix
T = eye( 4 );
T( 4, 1:3 ) = -t';

% Apply the translation
S = T * S * T';

% Return to original scale
S = quadric.normalize(S) .* abs(k);

end

