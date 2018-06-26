function [c, class] = matrixToCanonical(S)
% The canonical form of a quadric are the parameters:
%   cx, cy, cz, sx, sy, sz, alpha, beta, gamma, scale
%
% where [cx, cy, cz] if the center of the quadric, [sx, sy, sz] are the
% semi-axes, and [alpha, beta, gamma] are the angles in degrees.


% find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% translate to the center
Q = T * S * transpose(T);

% solve the eigenproblem
[evecs,evals] = svd(Q( 1:3, 1:3 ) / -Q( 4, 4 ));
radii = 1./diag(sqrt(evals));
sgns = sign( diag( evals ) );
radii = radii .* sgns;

% To rotate Q into canonical orientation:
%  evecs'*Q(1:3,1:3)*evecs

% Derive the angles
angles = rad2deg(rotm2eul(evecs));

% Derive the scale
scale = transpose(center)*S(1:3,1:3)*center - S(4,4);

% Assemble the canonical parameter vector
c = [center' radii' angles scale];

% Identify the class of the quadric
switch num2str(sign(diag(evals))','%d')
    case '111'
       class = 'ellipsoid';
    case {'11-1', '1-11', '-111'}
        class = 'hyperboloid of one sheet';
    case {'1-1-1', '-11-1', '-1-11'}
        class = 'hyperboloid of two sheets';
    otherwise
        class = 'degenerate case';
end

end

