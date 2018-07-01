syms A B C D E F G H I K
syms x y z
syms x0 y0 z0 x0d y0d z0d k0 n0
syms x1 y1 z1 x1d y1d z1d k1 n1

S = [A D E G; D B F H; E F C I; G H I K];
X = [x; y; z; 1];
X0 = [x0; y0; z0; 1];
X0d = [x0d; y0d; z0d; 1];
R0 = X0 + k1 * X0d;

X1 = [x1; y1; z1];
X1d = [x1d; y1d; z1d];
R1 = X1 + k1 * X1d;


% General equation for quadric surface:
eq2 = transpose(X) * S * X == 0;

% Ray Quadric intersection equation
eq11 = transpose([R0; 1]) * S * [R0; 1] == 0;

% Normal vector to a quadric
Q = 2 * S(1:3,1:4) * X1;
eq10 = Q / sqrt(sum(Q.^2));



% Snell's law
eq17 = cross(Q,n0*cross(X0d,Q)) == cross(Q,n1*cross(X1d,Q));
eq18 = solve(eq17,X1d);

X1d = [eq18.x1d; eq18.y1d; eq18.z1d];
