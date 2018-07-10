function S = normalize( S )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

scaleSign = sign(S(end,end));

% Adjust by the absolute value of K
S = S ./ abs(S(end));

% Restore the sign
S = S.*scaleSign;

end

