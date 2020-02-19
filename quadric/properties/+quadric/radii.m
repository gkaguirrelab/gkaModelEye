function r = radii(S)
% Return the semi-axes of the quadric, ordered in the x, y, z dimensions
%
% Syntax:
%  r = quadric.radii(S)
%
% Description:
%   A quadric has three semi-axis values. The semi-axis lengths are derived
%   from the matrix components A, B, C, and the scale component K. The
%   orientation of the semi-axis lengths to the x, y, z dimensions of the
%   space is influenced by the rotation of the quadric.
%
%   This routine returns the semi-axis lengths, ordered by their alignment
%   with the x, y, z dimensions.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   r                     - 3x1 vector providing semi-axes of the quadric
%                           that are most closely aligned with the [x y z]
%                           dimensions, respectively.
%
% Examples:
%{
    % Test that the routine returns the radii in the correct order
    p = perms([1 2 3]);
    for ii = 1:size(p,1)
        S = quadric.scale(quadric.unitSphere,p(ii,:));
        assert( isequal(quadric.radii(S),p(ii,:)'));
    end
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% Form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% Translate to the center
Q = T * S * transpose(T);

% Solve the eigenproblem
[~,evals] = svd(Q( 1:3, 1:3 ) / -Q( 4, 4 ));
r = 1./diag(sqrt(evals));
sgns = sign( diag( evals ) );
r = r .* sgns;

% The radius values are at this stage in canonical order (smallest to
% largest). Now re-order the valus so that the correspond to the actual
% x,y,z dimensions of the quadric.
dimensionRank = quadric.dimensionSizeRank(S);
r = r(dimensionRank);

end

