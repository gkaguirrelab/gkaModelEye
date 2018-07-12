function St = translate( S, Xt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%{
    S = quadric.scale(quadric.unitSphere,[5 4 3]);
    quadric.center(S)
    Xt = [-10; 2; -3];
    St = quadric.translate(S, Xt);
    quadric.center(St)
    Sprime = quadric.translate(St, -Xt);
    quadric.center(Sprime)

%}


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVars(S);

xt = Xt(1);
yt = Xt(2);
zt = Xt(3);

% Store the signs. These are needed to properly adjust the K element.
xs = sign(xt);
ys = sign(xt);
zs = sign(xt);

Gt = G - A*xt - D*yt - E*zt;
Ht = H - B*yt - D*xt - F*zt;
It = I - C*zt - E*xt - F*yt;
Kt = K + xs*A*xt^2 + ys*B*yt^2 + zs*C*zt^2 + E*xt*yt + D*xt*zt + F*yt*zt - G*xt - H*yt - I*zt;

St = quadric.varsToMatrix(A, B, C, D, E, F, Gt, Ht, It, Kt);

end

