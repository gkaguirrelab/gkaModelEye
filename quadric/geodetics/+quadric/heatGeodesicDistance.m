function [distance,xErrors] = heatGeodesicDistance(S,X0,X1,subdivisions)
% Find the geodesic distance between two points on a tri-axial ellipsoid
%
% Syntax:
%  [distance,startAngle,endAngle,geodeticPathCoords] = quadric.heatGeodesicDistance(S,G0,G1,X0,X1,subdivisions)
%
% Description:
%   Returns an estmiate of the geodesic distance between two points on an
%   un-rotated tri-axial ellipsoidal surface. In this approach, we begin
%   with a unit geodesic sphere created by subdividing a regular
%   icosahedron with normalised vertices. The position of the vertices of
%   the sphere are scaled to produce the specified ellipsoid. The heat
%   geodesic method is then used to estimate the distance between the
%   specified coordinates:
%
%       Crane, Keenan, Clarisse Weischedel, and Max Wardetzky. "Geodesics
%       in heat: A new approach to computing distance based on heat flow."
%       ACM Transactions on Graphics (TOG) 32.5 (2013): 1-11.
%
%   The heat geodesics code is taken from: 
%
%       http://www.numerical-tours.com/matlab/meshproc_7_geodesic_poisson/
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%   subdivisions          - Scalar. The resolution with which the
%                           ellipsoidal surface is sampled.
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
    % Numeric example provided by G. Panou in the original code
    S = quadric.scale(quadric.unitSphere,[0.015, 0.010, 0.009]);
    X0 = quadric.ellipsoidalGeoToCart([5; 5; 0],S);
    X1 = quadric.ellipsoidalGeoToCart([60; 120; 0],S);
    distance = quadric.heatGeodesicDistance(S,X0,X1);
    % Check the result against Panou's value
    assert( max(abs(distance - 0.0259)) < 1e-3 );
%}
%{
    % Distance from the fovea to the optic disc
    eye = modelEyeParameters('calcLandmarkFovea',true,'calcLandmarkOpticDisc',true);
    S = eye.retina.S;
    X0 = eye.landmarks.fovea.coords';
    X1 = eye.landmarks.opticDisc.coords';
    odf_distance = quadric.heatGeodesicDistance(S,X0,X1);
    outline = sprintf('Geodetic distance from the fovea to the optic disc: %2.2f mm\n',odf_distance);
    fprintf(outline);
%}


%% Handle arguments
% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

if nargin<4
    subdivisions=5;
end


%% Create the ellipsoid mesh

% Obtain the radii and center of the ellipsoid
radii = quadric.radii(quadric.alignAxes(S));
center = quadric.center( S );

% Create a unit sphere mesh, composed of vertices (X) and faces (F).
% Transpose these variables to match what the heat geodesics computations
% expect.
[X,F] = icosphere(subdivisions);
X = X';
F = F';

% Translate and scale the position of the vertices to produce our ellipsoid
X(1,:) = X(1,:) * radii(1) + center(1);
X(2,:) = X(2,:) * radii(2) + center(2);
X(3,:) = X(3,:) * radii(3) + center(3);

%% Add the start and end coordinates to the mesh
% This is currently commented out, as it didn't help much, so it wasn't
% worth adding this additional library function
% [ xErrors, ~, F, X] = ...
%     point2trimesh('Faces',F','Vertices',X','QueryPoints',[X0,X1]');
% X = X';
% F = F';

[~,startIdx] = min(vecnorm(X-X0));
[~,endIdx] = min(vecnorm(X-X1));


%% Calculate the distance
% The code that follows is taken directly from Gabriel Peyré's website,
% implementing the approach described by Crane and colleagues.
n = size(X,2);
m = size(F,2);

XF = @(i)X(:,F(i,:));
Na = cross( XF(2)-XF(1), XF(3)-XF(1) );

amplitude = @(X)sqrt( sum( X.^2 ) );
A = amplitude(Na)/2;

normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
N = normalize(Na);

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

GradMat = {};
for k=1:3
    GradMat{k} = dA*sparse(I,J,V(k,:),m,n);
end
Grad = @(u)[GradMat{1}*u, GradMat{2}*u, GradMat{3}*u]';

dAf = spdiags(2*A(:),0,m,m);
DivMat = {GradMat{1}'*dAf, GradMat{2}'*dAf, GradMat{3}'*dAf};
Div = @(q)DivMat{1}*q(1,:)' + DivMat{2}*q(2,:)' + DivMat{3}*q(3,:)';
Delta = DivMat{1}*GradMat{1} + DivMat{2}*GradMat{2} + DivMat{3}*GradMat{3};
cota = @(a,b)cot( acos( dot(normalize(a),normalize(b)) ) );
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


delta = zeros(n,1);
delta(startIdx) = 1;
t = 1e-4;
u = (Ac+t*DeltaCot)\delta;
g = Grad(u);
h = -normalize(g);
phi = Delta \ Div(h);
phi = phi - min(phi);

distance = phi(endIdx);

end



function [vv,ff] = icosphere(varargin)
%ICOSPHERE Generate icosphere.
% Create a unit geodesic sphere created by subdividing a regular
% icosahedron with normalised vertices.
%
%   [V,F] = ICOSPHERE(N) generates to matrices containing vertex and face
%   data so that patch('Faces',F,'Vertices',V) produces a unit icosphere
%   with N subdivisions.
%
%   FV = ICOSPHERE(N) generates an FV structure for using with patch.
%
%   ICOSPHERE(N) and just ICOSPHERE display the icosphere as a patch on the
%   current axes and does not return anything.
%
%   ICOSPHERE uses N = 3.
%
%   ICOSPHERE(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE.
%
%   Based on C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% Parse possible axes input
if nargin > 2
    error('Too many input variables, must be 0, 1 or 2.');
end
[cax,args,nargs] = axescheck(varargin{:});

n = 3; % default number of sub-divisions
if nargs > 0, n = args{1}; end % override based on input

% generate regular unit icosahedron (20 faced polyhedron)
[v,f] = icosahedron(); % size(v) = [12,3]; size(f) = [20,3];

% recursively subdivide triangle faces
for gen = 1:n
    f_ = zeros(size(f,1)*4,3);
    for i = 1:size(f,1) % for each triangle
        tri = f(i,:);
        % calculate mid points (add new points to v)
        [a,v] = getMidPoint(tri(1),tri(2),v);
        [b,v] = getMidPoint(tri(2),tri(3),v);
        [c,v] = getMidPoint(tri(3),tri(1),v);
        % generate new subdivision triangles
        nfc = [tri(1),a,c;
               tri(2),b,a;
               tri(3),c,b;
                    a,b,c];
        % replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_; % update 
end

% remove duplicate vertices
[v,b,ix] = unique(v,'rows'); clear b % b dummy / compatibility
% reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f),'rows');

switch(nargout)
    case 0 % no output
        cax = newplot(cax); % draw to given axis (or gca)
        showSphere(cax,f,v);
    case 1 % return fv structure for patch
        vv = struct('Vertices',v,'Faces',f,...
                    'VertexNormals',v,'FaceVertexCData',v(:,3));
    case 2 % return vertices and faces
        vv = v; ff = f;
    otherwise
        error('Too many output variables, must be 0, 1 or 2.');
end

end


function [i,v] = getMidPoint(t1,t2,v)
%GETMIDPOINT calculates point between two vertices
%   Calculate new vertex in sub-division and normalise to unit length
%   then find or add it to v and return index
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% get vertice positions
p1 = v(t1,:); p2 = v(t2,:);
% calculate mid point (on unit sphere)
pm = (p1 + p2) ./ 2;
pm = pm./norm(pm);
% add to vertices list, return index
i = size(v,1) + 1;
v = [v;pm];

end

function [v,f] = icosahedron()
%ICOSAHEDRON creates unit regular icosahedron
%   Returns 12 vertex and 20 face values.
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));

% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20
end
