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

% Normalize the Navarro values to have the standard form of -1 as a0
v = quadric.normalize(v);

% Convert from polynomial to matrix form
S = quadric.vecToMatrix(v);

% Rotate to canonical orientation
Srot = quadric.alignAxes( S );

% Translate to canonical center
[crot, class] = quadric.matrixToTransparent(Srot);
t = crot(7:9);
Srottrans = quadric.translate(S,t);


% Confirm that matrix can be recovered from the canonical form
Srecovered = quadric.transparentToMatrix(c);
assert(max(max(abs(Srecovered-S)))< 0.02);

% Examine the canonical form of the shifted quadric
cShift = quadric.matrixToTransparent(S);

% Obtain a function handle for the polynomial
F = quadric.returnPolynomialFunc(quadric.matrixToPolynomial(S));

% Plot the surface
gv = linspace(-15,15,100); % adjust for appropriate domain
[xx, yy, zz]=meshgrid(gv, gv, gv);
fig = figure;
isosurface(xx, yy, zz, F(xx, yy, zz), 0)
axis equal


foo=1;


