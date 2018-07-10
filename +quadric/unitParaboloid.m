function S = unitParaboloid
% opening along the first dimension

S = eye(4);
S(1,1)=0;
S(1,4)=-1;
S(4,1)=-1;
S(4,4)=0;

end