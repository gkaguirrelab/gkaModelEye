function S = rotate( S, r )
% Rotates a quadric surface by the specified Euler angles
%
% Syntax:
%  S = quadric.rotate( S, r )
%
% Description:
%   
%   S can be provided as a vector or matrix, and the returned quadric will
%   be of the same form.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   r                     - 1x3 vector of Euler angles
% Outputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface, having been subjected to rotation
%
% Examples:
%{
%}

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

