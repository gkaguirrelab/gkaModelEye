function N = surfaceNormal(S,X,side,surfaceTolerance)
% Returns a vector that is normal to a quadric surface at a coordinate
%
% Syntax:
%  N = surfaceNormal(S,X,side,surfaceTolerance)
%
% Description:
%   Given a quadric surface (S) and a coordinate on that surface, the
%   routine returns a vector that identifies the normal to the surface at
%   that coordiate.
%
% Inputs:
%   S                     - 4x4 quadratic surface matrix
%   X                     - 3x1 vector specifying a point on the surface
%   side                  - Scalar, taking the value of -1 or 1. The
%                           quadric surface has two sides at every point.
%                           This value determines if the "interior" or
%                           "exterior" surface normal is returned. For  an
%                           ellipsoid, setting side=-1 causes the routine
%                           to report the normal on the convex side, and
%                           setting side=+1 causes the routine to report
%                           the normal on the concave side. For a two-
%                           sheeted hyperboloid, with the ray traveling
%                           down the open axis, setting side = -1 returns
%                           the normal from the interior of one of the
%                           hyperboloid sheets, and thus the concave side.
%                           Default value is 1.
%   surfaceTolerance        Scalar. Defines the tolerance on the check to
%                           ensure that the passed X coordinate is on the
%                           quadric surface. Defaults to 1e-10 if not
%                           defined. If set to empty the check is skipped.                           
%
% Outputs:
%   N                     - 3x2 matrix that specifies the normal as a unit
%                           vector of the form [p; u], corresponding to
%                               R = p + t*u,
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t has an
%                           obligatory value of unity.
%
% Examples:
%{
    S = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    X = [0; 0.707106781186547; 0.707106781186547];
    N = quadric.surfaceNormal(S,X);
%}

% Keep the compiler happy by excluding prohibited calls
coder.extrinsic('warning')

% Handle incomplete input arguments
if nargin==2
    side=1;
    surfaceTolerance=1e-10;
end

if nargin==3
    surfaceTolerance=1e-10;
end

% Pre-allocate the output variables
N = nan(3,2);

% Clear the nan case
if any(isnan(X))
    return
end

% Check that the supplied point is on the surface of the quadric
if ~isempty(surfaceTolerance)
    funcS = quadric.vecToFunc(quadric.matrixToVec(S));
    % Decompose the coordinate
    x = X(1); y = X(2); z = X(3);
    if abs(funcS(x,y,z)) > surfaceTolerance
        warning('Passed coordinate is not on the quadric surface within tolerance (%f)',surfaceTolerance);
        return
    end
end

% Obtain the surface normal by taking the partial derivatives of Q with
% respect to x, y, and z. Performed here in matrix form.
Q = 2*S(1:3,:)*[X;1];
Q = Q/sqrt(dot(Q,Q));

% For a given point on the quadric surface, the normal could be expressed
% in either of two directions (i.e., directed towards the "exterior" of the
% surace or towards the "interior"). Whether the positive or negative value
% is controled by the side variable. If side==1, return -Q. In the case of
% an ellipsoid, this can be considered as reporting the normal directed
% towards the interior of the closed quadric, consistent with a ray having
% intersected the concave side of the surface.
if side==1
    Q=-Q;
end

% Express the surface normal as a unit vector arising from the point of
% intersection.
N(:,1)=X;
N(:,2)=Q;

end

