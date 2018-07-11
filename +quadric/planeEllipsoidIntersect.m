function [Aye,Bye,q1,q2,q3] = planeEllipsoidIntersect(S,p1,p2,p3)

n = cross(p1 - p2, p1 - p3);
d = p1(1)*n(1) + p1(2)*n(2) + p1(3)*n(3);
d = -d;


radii = quadric.radii(S);

[Aye,Bye,q1,q2,q3]=BektasFxn(n(1),n(2),n(3),d,3,4,5);

end




function[Aye,Bye,q1,q2,q3]=BektasFxn(A,B,C,D,a,b,c)

% This function computes the intersection an Ellipsoid and a Plane, 

% This source code was adapted from the implementation of the algorithms described in 
%   P. Klein, "On the Ellipsoid and Plane Intersection Equation," 
%   Applied Mathematics, Vol. 3 No. 11, 2012, pp. 1634-1640. doi: 10.4236/am.2012.311226.
 
% Plane equation          A.x + B.y + C.z + D = 0
% Ellipsoid equation    (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
% 

% Author: Sebahattin Bektas,Ondokuz Mayis University,2015 

%
% Inputs -----------------------------------------------------------------
%   A,B,C,D    :Plane equation's coefficients  
%                
%   a,b,c      : Semi-axes of ellipsoid 
%                
% Outputs ----------------------------------------------------------------
%   Aye:  Semi-major axis of intersection ellipse

%   Bye:  Semi-minor axis of intersection ellipse

%
% q1,q2,q3 : Cartesian coordinates of intersection ellipse's center


% This source code is given for free! However, I would be grateful if you refer 
% to corresponding article in any academic publication that uses 
% this code or part of it. Here are the corresponding references: 

%  BEKTAS, Sebahattin. Orthogonal distance from an ellipsoid. Bol. Ciênc. Geod. [online]. 2014, 
%  vol.20, n.4, pp. 970-983. ISSN 1982-2170.


%
car=A*B*C*D;
if car~=0
kx2=1/a^2 + (A^2)/(C^2*c^2);
ky2=1/b^2 + (B^2)/(C^2*c^2);
kxy=(2*A*B)/(C^2*c^2);
kx=+ (2*A*D)/(C^2*c^2);
ky=+ (2*B*D)/(C^2*c^2) ;
ksab=D^2/(C^2*c^2)- 1;
ParA=[kx2  kxy  ky2  kx  ky  ksab];
G=AtoG(ParA);

q1=G(1);q2=G(2);q3=(A*q1+B*q2+D)/-C;
    end

if C==0 & A~=0
ParA=[(1/b^2 + (B/A/a)^2) 0  1/c^2 (2*D*B/(A*a)^2)  0  (-1+(D/A/a)^2)];G=AtoG(ParA);

q2=G(1);q3=G(2);q1=(D+B*q2+C*q3)/-A;    
end

if C==0 & B~=0
ParA=[(1/a^2 + (A/B/b)^2) 0  1/c^2 (2*D*A/(B*b)^2)  0  (-1+(D/B/b)^2)];G=AtoG(ParA);

q1=G(1);q3=G(2);q2=(D+A*q1+C*q3)/-B;    
end

if A==0 & B==0 ,    q1=0;q2=0;q3=-D/C;end 
if D==0,    q1=0;q2=0;q3=0;

end 
 
quz=sqrt(q1^2+q2^2+q3^2);
n1=A/sqrt(A^2+B^2+C^2);
n2=B/sqrt(A^2+B^2+C^2);
n3=C/sqrt(A^2+B^2+C^2);
kap=(q1*A+B*q2+C*q3)/sqrt(A^2+B^2+C^2);
d= kap^2*(A^2+B^2+C^2)/(a^2*A^2+b^2*B^2+c^2*C^2);

ak=1;
bk=-(n1^2*(1/b^2+1/c^2)+n2^2*(1/c^2+1/a^2)+n3^2*(1/b^2+1/a^2));
ck=(n1/b/c)^2+(n2/a/c)^2+(n3/b/a)^2;
p=[ak bk ck];
kok=roots(p);
Bye=sqrt((1-d)/kok(1));
Aye=sqrt((1-d)/kok(2));

end

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