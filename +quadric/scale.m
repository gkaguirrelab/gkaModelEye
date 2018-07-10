function Sm = scale( S, s)

if length(s)==1
    s = [s;s;s];
end

Minv = eye(4);
Minv(1,1)=1/s(1);
Minv(2,2)=1/s(2);
Minv(3,3)=1/s(3);

Sm = Minv*S*Minv';

end

