function [distance,startAngle,endAngle,geodeticPathCoords] = poissonGeodesicDistance(S,G0,G1,X0,X1)
% Find the geodesic distance between two points on a tri-axial ellipsoid
%
% Syntax:
%  [distance,startAngle,endAngle,geodeticPathCoords] = poissonGeodesicDistance(S,X0,X1)
%
% Description:
%   Returns the geodesic distance between two points on the tri-axial
%   ellipsoidal surface. This is (effectively) the minimum length path on
%   the ellipsoidal surface that connects the two points. This approach
%   makes use of:
%
%       Crane, Keenan, Clarisse Weischedel, and Max Wardetzky. "Geodesics
%       in heat: A new approach to computing distance based on heat flow."
%       ACM Transactions on Graphics (TOG) 32.5 (2013): 152.
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   G0, G1                - 3x1 vectors that provide the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%
% Outputs:
%   distance              - Scalar. Distance of the geodetic between the
%                           two points.
%   startAngle, endAngle  - Scalars. The heading of the geodetic path, in
%                           degrees, relative to a line of constant omega
%                           in the ellipsoidal geodetic coordinate system
%                           on the ellipsoidal surface. 
%
% Examples:
%{
    eye = modelEyeParameters();
    S = eye.retina.S;
    G0 = [90;0;0];
    G1 = [25;100;0];
    [distance,startAngle] = quadric.poissonGeodesicDistance(S,G0,G1);
    tic
    [G1prime, distanceError, angleError] = quadric.geodesicByReckoning(S,G0,startAngle,distance);
    
[distance,startAngle,endAngle,geodeticPathCoords] = poissonGeodesicDistance(S,G0,G1,X0,X1)

%}

% resolution of the returned path coords
vertexDensity = 10;

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

S = quadric.alignAxes(S);

% If three input values were passed, convert the G0/G1 variables to X0/X1
% cartesian coordinates
if nargin==3
    % Obtain the ellipsoidal geodetic coordinates for the two points, and
    % convert to radians
    X0 = quadric.ellipsoidalGeoToCart( G0, S );
    X1 = quadric.ellipsoidalGeoToCart( G1, S );
end


% Obtain a set of points on the surface
betas = linspace(-90,0,vertexDensity);
omegas = linspace(-180,180,vertexDensity);
[BETAS,OMEGAS] = meshgrid(betas,omegas);
coordinates = [BETAS(:)'; OMEGAS(:)']';

% Obtain a Delauny triangulation for these coordinates
F = delaunay(coordinates);

% Obtain the Cartesian locations for the coordinates
X=[];
for ii=1:size(coordinates,1)
    X(ii,:) = quadric.ellipsoidalGeoToCart( coordinates(ii,:), S );
end

% Set X and F to have the expected row column order for the next process
X = X';
F = F';

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

fprintf('Should be 0: %e\n', norm(Delta-DeltaCot, 'fro')/norm(Delta, 'fro'));

% Set the index equal to 1 (the X0 point)
i = 25;
delta = zeros(n,1);
delta(i) = 1;

t = 1000;
u = (Ac+t*Delta)\delta;

options.face_vertex_color = u;
clf; plot_mesh(X,F,options);
axis('tight');
colormap parula(256);


% Compute the solution with explicit time stepping
t = .1;
u = (Ac+t*DeltaCot)\delta;

% Compute the gradient field.
g = Grad(u);

% Normalize it
h = -normalize(g);

% Integrate it back
phi = Delta \ Div(h);

figure
options.face_vertex_color = phi;
clf; plot_mesh(X,F,options);
axis('tight');
colormap parula(256);

end % panouGeodesicDistance


