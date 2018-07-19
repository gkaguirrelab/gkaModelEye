function Sm = scale( S, s)
% Scales the radii of a quardic surface by the passed scaling parameters
%
% Syntax:
%  Sm = quadric.scale( S, s)
%
% Description:
%   Returns a quadric surface that has had each radius scaled by the passed
%   set of scaling parameters. This may be used, for example, to transform
%   a sphere into an ellipsoid.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   s                     - 3x1 or 1x3 vector specifying the multiplicative
%                           scaling to be applied to each semi-radius, in
%                           order. If passed as a scalar, then the single
%                           value is applied to all three dimensions.
%
% Outputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Examples:
%{
    % Create an ellipsoid
    S = quadric.unitSphere();
    S = quadric.scale(S,[3 4 5]);
    % Report the radii
    quadric.radii(S)
%}


% If the quadric surface was passed in vector form, convert to matrix
returnVecFlag = false;
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
    returnVecFlag = true;
end

% If s was passed as a scalar, expand to three dimensions
if length(s)==1
    s = [s;s;s];
end

% Apply the scaling
Minv = eye(4);
Minv(1,1)=1/s(1);
Minv(2,2)=1/s(2);
Minv(3,3)=1/s(3);
Sm = Minv*S*Minv';

% Return a vector if that was the original input
if returnVecFlag
    Sm = quadric.matrixToVec(Sm);
end


end

