function St = translate( S, Xt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVariables(S);

xt = Xt(1);
yt = Xt(2);
zt = Xt(3);

Gt = G - A*xt - D*yt - E*zt;
Ht = H - B*yt - D*xt - F*zt;
It = I - C*zt - E*xt - F*yt;
Kt = K + A*xt^2 + B*yt^2 + C*zt^2 + E*xt*yt + D*xt*zt + F*yt*zt - G*xt - H*yt - I*zt;

St = quadric.variablesToMatrix(A, B, C, D, E, F, Gt, Ht, It, Kt);

end

