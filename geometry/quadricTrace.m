syms A B C D E F G H I K
syms x y z
syms x0 y0 z0 x0d y0d z0d k0 n0
syms x1 y1 z1 x1d y1d z1d k1 n1

S = [A D F G; D B E H; F E C I; G H I K];
X = [x; y; z; 1];
X0 = [x0; y0; z0];
X0d = [x0d; y0d; z0d];
R0 = X0 + k1 * X0d;

X1 = [x1; y1; z1];
X1d = [x1d; y1d; z1d];
R1 = X1 + k1 * X1d;


Q = 2 * S(1:3,1:4) * X;
Qn = Q / sqrt(sum(Q.^2));

% General equation for quadric surface:
eq2 = transpose(X) * S * X == 0;

% Ray Quadric intersection equation
eq11 = transpose([R0; 1]) * S * [R0; 1] == 0;

% Snell's law
eq17 = cross(Q,n0*cross(X0d,Q)) == cross(Q,n1*cross(X1d,Q));
eq18 = solve(eq17,X1d);

X1d = [eq18.x1d; eq18.y1d; eq18.z1d];

% Translation of a quadric surface
syms xt yt zt
Xt = [xt; yt; zt; 1];
transpose(X-Xt)*S*(X-Xt);



function S = surfaceMatrix(p)

% The implicit form of a second-order (quadric) surface:
%   S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
%               2Dxy + 2Exy + 2Fxz +
%               2Gx + 2Hy + 2Iz + K == 0
%

A = p(1);
B = p(2);
C = p(3);
D = p(4);
E = p(5);
F = p(6);
G = p(7);
H = p(8);
I = p(9);
K = p(10);

S = [A D F G; D B E H; F E C I; G H I K];

end
