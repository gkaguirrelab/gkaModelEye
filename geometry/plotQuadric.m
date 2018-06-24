% Plot the canonical Navarro corneal surface


% Navarro 2006, Table 1, mean cornea ellipsoid parameters
a11 = 1; a22 = 1.0318; a33 = 0.536; a12 = -0.0006; a13 = -0.038;
a23 = -0.0147; a1 = -.06; a2 = -0.115; a3 = 13.546; a0 = -23.25;

% Place the Navarro parameters into a polynomial vector. Navarro defines his
% parameters w.r.t. the polynomial equation as:
%
%   0 =
%     a11x^2 + a22y^2 + a33z^2 + a12xy + a13xz + a23yz + a1x + a2y +a3z+a0
%
% To convert to the standard notation, need to account for
% the factor of 2. Note the cross-terms are in the order [xy, xz, yz].
v = [a11 a22 a33 a12/2 a13/2 a23/2 a1/2 a2/2 a3/2 a0];
   
% Convert from polynomial to matrix form
S = quadricPolynomialToMatrix(v);

% Convert the quadric from matrix to canonical form
c = quadricMatrixToCanonical(S);

% Convert the canonical form back to matrix
Srecovered = quadricCanonicalToMatrix(c);
vrecovered = quadricMatrixToPolynomial(Srecovered);

% Confirm that we can transform these representations accurately
assert(max(max(S-Srecovered)) < 1e-10);
assert(max(max(v-vrecovered)) < 1e-10);


% Obtain a function handle for the polynomial
F = quadricPolynomialFunc(v);

% Plot the surface
gv = linspace(-4,4,100); % adjust for appropriate domain
[xx, yy, zz]=meshgrid(gv, gv, gv);
fig = figure;
isosurface(xx, yy, zz, F(xx, yy, zz), 0)



foo=1;


function S = quadricCanonicalToMatrix(c)


center = c(1:3);
radii = c(4:6);
angles = c(7:9);
scale = c(10);

eulerAngles = deg2rad(angles.*([1 -1 1])-[90 0 -180]);

D = diag(radii);
Ainv = eul2rotm(eulerAngles)*D;
A = inv(Ainv);
Q = A'*A*scale;
Q(4,4) = -scale;

% Translation
T = eye( 4 );
T( 4, 1:3 ) = center;
sgns = ones(4,4);
sgns(4,1:3) = sign(center);
sgns(1:3,4) = sign(center);

S = T * Q * T' .* sgns;

end


function [c, class] = quadricMatrixToCanonical(S)
% The canonical form of a quadric are the parameters:
%   cx, cy, cz, sx, sy, sz, alpha, beta, gamma, scale
%
% where [cx, cy, cz] if the center of the quadric, [sx, sy, sz] are the
% semi-axes, and [alpha, beta, gamma] are the angles in degrees.


% find the center of the ellipsoid
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

% Derive the angles
eulerAngles = rotm2eul(evecs);

% Conver the angles into the style used by Navarro
angles = (rad2deg(eulerAngles)+[90 0 -180]).*([1 -1 1]);

% Derive the scale
scale = transpose(center)*S(1:3,1:3)*center - S(4,4);

c = [center' radii' angles scale];

end


function F = quadricPolynomialFunc(v)

% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%

F = @(xx,yy,zz) v(1)*xx.^2 + v(2)*yy.^2 + v(3)*zz.^2 + 2*v(4)*xx.*yy + 2*v(5)*xx.*zz + 2*v(6)*yy.*zz + 2*v(7)*xx + 2*v(8)*yy + 2*v(9)*zz + v(10);

end


function v = quadricMatrixToPolynomial(S)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%

v(1) = S(1,1);
v(2) = S(2,2);
v(3) = S(3,3);
v(4) = S(1,2);
v(5) = S(1,3);
v(6) = S(2,3);
v(7) = S(1,4);
v(8) = S(2,4);
v(9) = S(3,4);
v(10) = S(4,4);

end


function S = quadricPolynomialToMatrix(v)
% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exz + 2Fyz +
%               2Gx + 2Hy + 2Iz + K == 0
%
% Note that the order of the cross-terms is xy, xz, yz
%

A = v(1);
B = v(2);
C = v(3);
D = v(4);
E = v(5);
F = v(6);
G = v(7);
H = v(8);
I = v(9);
K = v(10);

S = [A D E G; D B F H; E F C I; G H I K];

end
