function S = translate( S, Xshift )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

Xt = [Xshift; 1];
S = Xt'.*S.*Xt;

end

