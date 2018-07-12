function St = translate( S, Xt )
%{
    % Define a quadric and a translation vector
    S = quadric.scale(quadric.unitSphere,[5 4 3]);
    Xt = [-10; 0; 0];

    % Translate the quadric out and back
    Sprime = quadric.translate(quadric.translate(S, Xt), -Xt);
    
    % Confirm recovery of original values
    assert(max(max(abs(S-Sprime))) < 1e-20);
%}

returnVecFlag = false;

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
    returnVecFlag = true;
end

% Obtain the variable form
[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVars(S);

% decompose the translation vector
xt = Xt(1);
yt = Xt(2);
zt = Xt(3);

% Adjust the terms
Gt = G - A*xt - D*yt - E*zt;
Ht = H - B*yt - D*xt - F*zt;
It = I - C*zt - E*xt - F*yt;
Kt = K + A*xt^2 + B*yt^2 + C*zt^2 + 2*D*xt*yt + 2*E*xt*zt + 2*F*yt*zt - 2*G*xt - 2*H*yt - 2*I*zt;

% Assemble the translated quadric matrix
St = quadric.varsToMatrix(A, B, C, D, E, F, Gt, Ht, It, Kt);

% Return a vector if that was the original input
if returnVecFlag
    St = quadric.MatrixToVec(St);
end

end

