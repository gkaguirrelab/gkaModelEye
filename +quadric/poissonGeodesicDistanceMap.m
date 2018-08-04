function [phi, X] = poissonGeodesicDistanceMap(S,G0,X0, nVertices)
% Compute a map of geodesic distance on a tri-axial ellipsoid
%
% Syntax:
%  [distance,startAngle,endAngle,geodeticPathCoords] = poissonGeodesicDistance(S,X0,X1)
%
% Description:
%   Returns the geodesic distance between two points on the tri-axial
%   ellipsoidal surface. This is (effectively) the minimum length path on
%   the ellipsoidal surface that connects the two points. This routine
%   is based on the approach of Crane and colleagues:
%
%       Crane, Keenan, Clarisse Weischedel, and Max Wardetzky. "Geodesics
%       in heat: A new approach to computing distance based on heat flow."
%       ACM Transactions on Graphics (TOG) 32.5 (2013): 152.
%
%   With the MATLAB implementation taken from:
%
%       G. Peyré, The Numerical Tours of Signal Processing - Advanced
%       Computational Signal and Image Processing IEEE Computing in Science
%       and Engineering, vol. 13(4), pp. 94-97, 2011.
%       http://www.numerical-tours.com/matlab/meshproc_7_geodesic_poisson/
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   G0                    - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0                    - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   nVertices             - Scalar. The requested number of vertices used
%                           to triangulate the ellipsoidal surface. The
%                           precise number of vertices obtained will differ
%                           somewhat from this value.
%
% Outputs:
%   phi                   - 1xnVertices vector. Distance of the geodetic
%                           between X0 and each of the points listed in X. 
%   X                     - 3xnVertices matrix. The set of points in 
%                           Cartesian coordinates that triangulate the
%                           ellipsoidal surface.
%
% Examples:
%{
    eye = modelEyeParameters();
    S = eye.retina.S;
    X0 = eye.axes.visual.coords;
    [phi, X] = quadric.poissonGeodesicDistanceMap(S,[],X0 );
    c = jet();
    nColors = size(c,1);
    figure
    hold on
    for ii=1:size(X,2)
        colorTriple = c(round((phi(ii)./max(phi))*(nColors-1)+1),:);
    	plot3(X(1,ii),X(2,ii),X(3,ii),'.','MarkerSize',20,'Color',colorTriple);
    end
%}


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

S = quadric.alignAxes(S);

% If three input values were passed, assume that the X0/X1 variables are
% empty, and compute these from the passed G0/G1 geodetic coordinates
if nargin==2
    % Obtain the ellipsoidal geodetic coordinates for the two points, and
    % convert to radians
    X0 = quadric.ellipsoidalGeoToCart( G0, S );
    % Set the nVertices
    nVertices = 10000;
end

if nargin==3
    % Set the nVertices
    nVertices = 10000;
end

% Create a triangulation of a unit sphere with nearly equal areas. This step makes use of the spheretri
% respository, written by Peter Gagarinov, PhD
%   https://github.com/pgagarinov/spheretri
[vMat, fMat] = spheretri(nVertices);

% Project the vertices from the unit sphere to the ellipsoid
unitSphere = quadric.unitSphere;
X=[];
for ii = 1: size(vMat,1)
    geoCoord = quadric.cartToParametricGeo(vMat(ii,:)',unitSphere);
    X(ii,:) = quadric.parametricGeoToCart(geoCoord',S);
end

% Find the index of the point in the triangulation that is closest to X0
[~,idx] = min(sum((X-X0).^2,2));

% Set X and F to have the expected row column order for the next process
X = X';
F = fMat';


%% The following code is taken directly from Numerical Tours

% Number n vertices and number m faces
n = size(X,2);
m = size(F,2);

% Callback to get the coordinates of all the vertex of index i=1,2,3 in all faces
XF = @(i)X(:,F(i,:));

% Compute un-normalized normal
Na = cross( XF(2)-XF(1), XF(3)-XF(1) );

% Compute the area of each face as half the norm of the cross product.
amplitude = @(X)sqrt( sum( X.^2 ) );
A = amplitude(Na)/2;

% Compute the set of unit-norm normals to each face.
normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
N = normalize(Na);

% Populate the sparse entries of the matrices
I = []; J = []; V = []; % indexes to build the sparse matrices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % vector N_f^e_i
    wi = cross(XF(t)-XF(s),N);
    % update the index listing
    I = [I, 1:m];
    J = [J, F(i,:)];
    V = [V, wi];
end
dA = spdiags(1./(2*A(:)),0,m,m);

% Compute gradient
GradMat = {};
for k=1:3
    GradMat{k} = dA*sparse(I,J,V(k,:),m,n);
end

% Gradient operator
Grad = @(u)[GradMat{1}*u, GradMat{2}*u, GradMat{3}*u]';

% Compute divergence matrices as transposed of grad for the face area inner product.
dAf = spdiags(2*A(:),0,m,m);
DivMat = {GradMat{1}'*dAf, GradMat{2}'*dAf, GradMat{3}'*dAf};

% Div operator
Div = @(q)DivMat{1}*q(1,:)' + DivMat{2}*q(2,:)' + DivMat{3}*q(3,:)';

% Laplacian operator as the composition of grad and div.
Delta = DivMat{1}*GradMat{1} + DivMat{2}*GradMat{2} + DivMat{3}*GradMat{3};

% Cotan of an angle between two vectors.
cota = @(a,b)cot( acos( dot(normalize(a),normalize(b)) ) );

% Compute cotan weights Laplacian.
I = []; J = []; V = []; % indexes to build the sparse matrices
Ia = []; Va = []; % area of vertices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % adjacent edge
    ctheta = cota(XF(s)-XF(i), XF(t)-XF(i));
    % ctheta = max(ctheta, 1e-2); % avoid degeneracy
    % update the index listing
    I = [I, F(s,:), F(t,:)];
    J = [J, F(t,:), F(s,:)];
    V = [V, ctheta, ctheta];
    % update the diagonal with area of face around vertices
    Ia = [Ia, F(i,:)];
    Va = [Va, A];
end
% Aread diagonal matrix
Ac = sparse(Ia,Ia,Va,n,n);
% Cotan weights
Wc = sparse(I,J,V,n,n);
% Laplacian with cotan weights.
DeltaCot = spdiags(full(sum(Wc))', 0, n,n) - Wc;

% Set the index equal to the X0 point
delta = zeros(n,1);
delta(idx) = 1;

% Compute the solution with explicit time stepping
t = .1;
u = (Ac+t*DeltaCot)\delta;

% Compute the gradient field.
g = Grad(u);

% Normalize it
h = -normalize(g);

% Integrate it back
phi = Delta \ Div(h);

% Set the distance at the starting point to zero, and transpose to match
% the orientation of X.
phi = (phi - phi(idx))';


end % panouGeodesicDistance


