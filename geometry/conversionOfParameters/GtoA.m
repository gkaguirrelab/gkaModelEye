function ParA = GtoA(ParG);

%  Conversion of geometric parameters of an ellipse
%     ParG = [Center(1:2), Axes(1:2), Angle]'
%  to algebraic parameters
%  ParA = [A,B,C,D,E,F]'-  parameter vector of the fitting
%  conic:  Ax^2 + Bxy + Cy^2 +Dx + Ey + F = 0
%
%  Copyright 2011 Hui Ma

c = cos(ParG(5));  s = sin(ParG(5));
a = ParG(3);  b = ParG(4);
Xc = ParG(1);  Yc = ParG(2);
P = (c/a)^2 + (s/b)^2;
Q = (s/a)^2 + (c/b)^2;
R = 2*c*s*(1/a^2 - 1/b^2);
ParA = [P; R; Q; -2*P*Xc-R*Yc; -2*Q*Yc-R*Xc; P*Xc^2+Q*Yc^2+R*Xc*Yc-1];

end

