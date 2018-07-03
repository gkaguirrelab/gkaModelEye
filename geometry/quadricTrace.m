syms A B C D E F G H I K
syms x y z
syms px py pz ux uy uz t

S = [A D E G; D B F H; E F C I; G H I K];
X = [x; y; z];
p = [px; py; pz];
u = [ux; uy; uz];
R = p + t * u;


% General equation for quadric surface:
eq2 = transpose([X;1]) * S * [X;1] == 0;

% Ray Quadric intersection equation
eq11 = transpose([R0; 1]) * S * [R0; 1] == 0;

% Normal vector to a quadric
Q = 2 * S(1:3,1:4) * X1;
eq10 = Q / sqrt(sum(Q.^2));



% Snell's law
eq17 = cross(Q,n0*cross(X0d,Q)) == cross(Q,n1*cross(X1d,Q));
eq18 = solve(eq17,X1d);

X1d = [eq18.x1d; eq18.y1d; eq18.z1d];
