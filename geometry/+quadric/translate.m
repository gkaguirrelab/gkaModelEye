function S = translate( S, t )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% translate to the center
Q = T * S * transpose(T);


Xt = [t; 1];
S = Xt'.*S.*Xt;

end

