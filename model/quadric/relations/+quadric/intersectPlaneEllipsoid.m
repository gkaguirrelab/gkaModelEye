function [Aye,Bye,q1,q2,q3] = intersectPlaneEllipsoid(A,B,C,D,a,b,c)
% Ellipse parameters for a plane-ellipsoid intersection
%
% Syntax:
%  [Aye,Bye,q1,q2,q3] = quadric.intersectPlaneEllipsoid(A,B,C,D,a,b,c)
%
% Description:
%
%   Author: Sebahattin Bektas,Ondokuz Mayis University,2015
%
%   This function computes the intersection an ellipsoid and a plane, based
%   upon the algorithms described in
%
%       P. Klein, "On the Ellipsoid and Plane Intersection Equation,"
%       Applied Mathematics, Vol. 3 No. 11, 2012, pp. 1634-1640. doi:
%       10.4236/am.2012.311226.
%
%	Plane equation          A.x + B.y + C.z + D = 0
%	Ellipsoid equation    (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
%
% Inputs:
%	A,B,C,D               - Scalar. Coeffeicients of the plane equation
%   a,b,c                 - Scalar. Semi-axes of ellipsoid
%
% Outputs:
%   Aye                   - Scalar. Semi-major axis of intersection ellipse
%   Bye                   - Scalar. Semi-minor axis of intersection ellipse
%   q1,q2,q3              - Scalr. Cartesian coordinates of center of the
%                           ellipse
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


% AtoG
% Convert conic section equation from ABCDEF form to {center, axes, angle}.
% Minor modifications made to the original AtoG posted by Hui Ma
%  Note: Ma claimed copyright, but that is automatically released to the
% public upon his posting to FileExchange
%
%  ParA = [A,B,C,D,E,F]'-  parameter vector of the conic:  Ax^2 + Bxy + Cy^2 +Dx + Ey + F = 0
%  to geometric parameters  ParG = [Center(1:2), Axes(1:2), Angle]'
%
% The Angle value is in radians
%  code is:  1 - ellipse
%            2 - hyperbola
%            3 - parabola
%           -1 - degenerate cases
%            0 - imaginary ellipse
%            4 - imaginary parelle lines
%
%
function [ParG,code] = AtoG(ParA)
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
% some people never learn...
Angle = mod(Angle,pi);
% while Angle > pi
%     Angle = Angle - pi;
% end
% while Angle < 0
%     Angle = Angle + pi;
% end
ParG = [XYcenter; Axes; Angle];
end
