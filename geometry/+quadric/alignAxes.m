function S = alignAxes(S)

[evecs,~] = svd(S( 1:3, 1:3 ) / -S( 4, 4 ));
S(1:3,1:3) = evecs'*S(1:3,1:3)*evecs;

end