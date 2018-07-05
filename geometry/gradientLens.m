b = 4.8;
c0 = 1.418;
c1 = 1.371 - c0;
dL1 = 1.44;
dL2 = 2.16;

% lens front
N0params = [c0+c1, -2*c1/dL1, c1/(dL1^2)];
N1params = [c1/b^2];
N0 = @(Z) N0params(1) + N0params(2)*Z + N0params(3)*Z^2;
N1 = @(Z) N1params(1);

nFrontAnon = @(X,Y,Z) N0(Z)+N1(Z)*(X^2+Y^2);

% lens back
N0params = [c0, 0, c1/(dL2^2)];
N1params = [c1/b^2];
N0 = @(Z) N0params(1) + N0params(2)*Z + N0params(3)*Z^2;
N1 = @(Z) N1params(1);

nBackAnon = @(X,Y,Z) N0(Z)+N1(Z)*(X^2+Y^2);

nFront = func2str(nFrontAnon);
nBack = func2str(nBackAnon);

nfront = @(X,Y,Z) 1.371000000000000 + 0.065277777777778*Z  -0.022665895061728*Z^2 -0.002039930555556*(X^2+Y^2);
nBack = @(X,Y,Z) 1.418000000000000 + 0*Z  -0.010073731138546*Z^2 -0.002039930555556*(X^2+Y^2);

nFrontQuadric = [-0.022665895061728 0 0 0.065277777777778; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0.065277777777778 0 0 1.371000000000000];
nBackQuadric = [-0.010073731138546 0 0 0.0; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0 0 0 1.418000000000000];
