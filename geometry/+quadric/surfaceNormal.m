function N = surfaceNormal(S,X)
%
% Inputs:
%   S                     - 4x4 quadratic surface matrix
%   X                     - 3x1 vector specifying a point on the surface
%
% Outputs:
%   R                     - 3x2 matrix that specifies the normal as a unit
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u,
%                           where p is vector origin, d is the direction
%                           expressed as a unit steo, and t has an
%                           obligatory value of unity.
%
% Examples:
%{
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    X = [0; 0.707106781186547; 0.707106781186547];
    N = quadric.surfaceNormal(S,X);
%}

% Pre-allocate the output variables
R = nan(3,2);

% Decompose the coordinate
x = X(1); y = X(2); z = X(3);

% Test that the supplied point is on the surface of the quadric
funcS = quadric.vecToFunc(quadric.matrixToVec(S));
if abs(funcS(x,y,z)) > 1e-10
    warning('Passed coordinate is not on the surface of the passed quadric');
    return
end

% Obtain the surface normal by taking the partial derivatives of Q with
% respect to x, y, and z. Performed here in matrix form.
Q = 2*S(1:3,:)*[X;1];
Q = Q/sqrt(dot(Q,Q));

% Express the surface normal as a unit vector arising from the point of
% intersection.
N(:,1)=X;
N(:,2)=Q;

end

