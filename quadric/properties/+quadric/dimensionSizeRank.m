function dimensionRank = dimensionSizeRank(S)
% Returns the ordering of dimensions by the largest to smallest semi-axes
%
% Syntax:
%  dimensionRank = quadric.dimensionSizeRank(S)
%
% Description:
%   A quadric has three semi-axis values (A B C), defined in three
%   dimensions (x, y, z). For any given quadric, which may further be
%   rotated, the axis values [A B C] will be most closely oriented with a
%   particular dimension. This routine returns a list of the dimensions,
%   ranked (largest to smallest) by the size of the semi-axis that is
%   closest to it. For example, if a given quadric has it's longest
%   semi-axis most closely oriented with the y dimension, and it's shortest
%   semi-axis most closely oriented wth the z dimension, then the routine
%   returns the dimension rank [2 1 3].
%
%   As a convention, the radius along the (e.g.) x dimension is the radius
%   that correspond to the axis of the quadric surface that has an angle
%   between -45 and 45 degrees w.r.t. to the x axis.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   dimensionOrder        - 1x3 vector providing the rank ordering of the
%                           space determined by the lengths of the
%                           semi-axes of the quadric.
%
% Examples:
%{
    S = quadric.scale(quadric.unitSphere,[5 3 4]);
    quadric.dimensionSizeRank(S)
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Solve the eigenproblem
[evecs,~] = svd(-S( 1:3, 1:3 ) );

% Obtain the Euler angles
EulerAngles = rad2deg(rotm2eul(evecs));

% Identify the orientation of the axes
thisOrientation = mod(fix((abs(EulerAngles)+45)./90),2);

% Define a mapping of orientations to radii order
orientations = {[0 0 0],[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[1 1 1]};
orders = {[1 2 3],[1 3 2],[3 2 1],[2 1 3],[2 3 1],[3 1 2],[3 2 1]};

% Return the rank ordering of the dimensions
dimensionRank = orders{cellfun(@(x) isequal(x,thisOrientation),orientations)};

end