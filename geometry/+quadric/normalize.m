function S = normalize( S )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


% Adjust by the absolute value of K
S = S ./ abs(S(end,end));

end

