function [ParG,code] = AtoG(ParA)

%  Conversion of algebraic parameters 
%  ParA = [A,B,C,D,E,F]'-  parameter vector of the conic:  Ax^2 + Bxy + Cy^2 +Dx + Ey + F = 0
%  to geometric parameters  ParG = [Center(1:2), Axes(1:2), Angle]'

%  code is:  1 - ellipse  
%            2 - hyperbola 
%            3 - parabola
%           -1 - degenerate cases  
%            0 - imaginary ellipse  
%            4 - imaginary parelle lines
%
%
%  Copyright 2011 Hui Ma

tolerance1 = 1.e-10;
tolerance2 = 1.e-20;

% ParA = ParA/norm(ParA);
if (abs(ParA(1)-ParA(3)) > tolerance1)
    Angle = atan(ParA(2)/(ParA(1)-ParA(3)))/2;
else
    Angle = pi/4;
end

c = cos(Angle);  s = sin(Angle);
Q = [c s; -s c];
M = [ParA(1)  ParA(2)/2;  ParA(2)/2  ParA(3)];

D = Q*M*Q';
N = Q*[ParA(4); ParA(5)];
O = ParA(6);

if (D(1,1) < 0) && (D(2,2) < 0)
    D = -D;
    N = -N;
    O = -O;
end

UVcenter = [-N(1,1)/2/D(1,1); -N(2,1)/2/D(2,2)];

free = O - UVcenter(1,1)*UVcenter(1,1)*D(1,1) - UVcenter(2,1)*UVcenter(2,1)*D(2,2);

% if the determinant of [A B/2 D/2;B/2 C E/2;D/2 E/2 F]is zero 
% and if K>0,then it is an empty set;
% otherwise the conic is degenerate

Deg =[ParA(1),ParA(2)/2,ParA(4)/2;...
     ParA(2)/2,ParA(3),ParA(5)/2;...
     ParA(4)/2,ParA(5)/2,ParA(6)];
K1=[ParA(1),ParA(4)/2;ParA(4)/2 ParA(6)];
K2=[ParA(3),ParA(5)/2;ParA(5)/2 ParA(6)];
K = det(K1)+ det(K2);

if (abs(det(Deg)) < tolerance2)
    if (abs(det(M))<tolerance2) &&(K > tolerance2)
        code = 4;  % empty set(imaginary parellel lines)
    else
        code = -1; % degenerate cases
    end
else
    if (D(1,1)*D(2,2) > tolerance1)
        if (free < 0)
            code = 1; % ellipse
        else
            code = 0; % empty set(imaginary ellipse)
        end
    elseif (D(1,1)*D(2,2) < - tolerance1)
        code = 2;  % hyperbola
    else
        code = 3;  % parabola
    end
end

XYcenter = Q'*UVcenter;
Axes = [sqrt(abs(free/D(1,1))); sqrt(abs(free/D(2,2)))];

if code == 1 && Axes(1)<Axes(2)
    AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
    Angle = Angle + pi/2;
end

if code == 2 && free*D(1,1)>0
    AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
    Angle = Angle + pi/2;
end

while Angle > pi
    Angle = Angle - pi;
end
while Angle < 0
    Angle = Angle + pi;
end

ParG = [XYcenter; Axes; Angle];

end
