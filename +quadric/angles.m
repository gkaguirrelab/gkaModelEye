function angles = angles(S)
% Returns the angles of a quadric surface with respect to the axes
%
% Syntax:
%  angles = quadric.angles(S)
%
% Description:
%   Still figuring out the angles. This is an approach that yields angles
%   similar to that reported by Navarro for the average normal corneal
%   surface w.r.t. the keratometric axis.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   angles                - 1x3 vector containing the angles of rotation in
%                           degrees.
%   
% Examples:
%{
    % Obtain the rotation angles of Navarro JOSA A 23.2 (2006) cornea
    a11 = 1; a22 = 1.0318; a33 = 0.536; a12 = -0.0006; a13 = -0.038;
    a23 = -0.0147; a1 = -.06; a2 = -0.115; a3 = 13.546; a0 = -23.25;
    S = quadric.vecToMatrix([a11 a22 a33 a12/2 a13/2 a23/2 a1/2 a2/2 a3/2 a0]);
    % Confirm the semi-radii lengths match the size and order of Navarro
    r = quadric.radii(S);
    assert(max(abs(r-[10.43; 10.27; 14.26]))<0.01)
    % Confirm that the angles alpha, beta, gamma, match Navarro
    % [0.85, -2.35, 0.02]
    a = quadric.angles(S);
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Solve the eigenproblem
[evecs,~] = svd(-S( 1:3, 1:3 ) );

% Derive the angles
alfax=atan2(-evecs(3,2),evecs(3,3));
alfay=acos(evecs(3,3)/cos(alfax));
alfaz=atan2(-evecs(2,1),evecs(1,1));
angles=180/pi.*[alfax,alfay,alfaz];

end

