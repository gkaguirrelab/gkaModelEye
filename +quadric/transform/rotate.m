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
%   r                     - 1x3 vector of Euler angles in degrees
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

% Find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% Form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% Translate to the center
Q = T * S * transpose(T);

% Convert deg to rad
ro = pi/180;

% Create the rotation matrix
Rx = [1 0 0; 0 cos(r(1)*ro) sin(r(1)*ro); 0 -sin(r(1)*ro) cos(r(1)*ro)];
Ry = [cos(r(2)*ro) 0 -sin(r(2)*ro); 0 1 0; sin(r(2)*ro) 0 cos(r(2)*ro)];
Rz = [cos(r(3)*ro) sin(r(3)*ro) 0; -sin(r(3)*ro) cos(r(3)*ro) 0; 0 0 1];
Rmat = Rx*Ry*Rz;

% Apply the rotation rotation
Q(1:3,1:3) = Rmat'*Q(1:3,1:3)*Rmat;

% Restore the translation
T = eye( 4 );
T( 4, 1:3 ) = -center';
S = T * Q * transpose(T);

% Restore the original scale.
if abs(Sscale) > realmin
    S = (S./S(4,4)).*Sscale;
end

% Return a vector if that was the original input
if returnVecFlag
    S = quadric.matrixToVec(S);
end

end

