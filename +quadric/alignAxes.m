function [S, evecs] = alignAxes(S)
% Rotates a quadric surface so that it is in canonical form
%
% Syntax:
%  [S, evecs] = quadric.alignAxes(S)
%
% Description:
%   Given a quadric surface S, this routine will rotate the surface so that
%   the semi-axes are aligned with the [x y z] dimensions of the Cartesian
%   space, with the smallest semi-axis aligned with the x dimension and the
%   largest semi-axis aligned with the z dimension.
%
%   S can be provided as a vector or matrix, and the returned quadric will
%   be of the same form.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
% Outputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface, now in canonical form
%   evecs                 - The rotation matrix used to place the quadric
%                           in canonical form
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

% Here and elsewhere, detect the case in which the scale is close to zero,
% in which case no adjustment is made. This can occur when the quadric is
% translated so that an apex is positioned at the origin.
if Sscale < realmin 
    Sscale = 1;
end

% Find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% Form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% Translate to the center
Q = T * S * transpose(T);

% Store the translated scale
Qscale = Q(4,4);

% Here and elsewhere, detect the case in which the scale is close to zero,
% in which case no adjustment is made. This can occur when the quadric is
% translated so that an apex is positioned at the origin.
if Qscale < realmin 
    Qscale = 1;
end

% Solve the eigenproblem.
[evecs,~] = svd(Q(1:3,1:3) / -Qscale);

% Apply the inverse rotation
Q(1:3,1:3) = evecs'*Q(1:3,1:3)*evecs;

% Restore the Q scale
Q = (Q./Q(4,4)).*Qscale;

% Restore the translation
T = eye( 4 );
T( 4, 1:3 ) = -center';
S = T * Q * transpose(T);

% Restore the original scale.
if Sscale > 1e-12
    S = (S./S(4,4)).*Sscale;
end

% Return a vector if that was the original input
if returnVecFlag
    S = quadric.matrixToVec(S);
end


end