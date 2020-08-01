function X = mostAnteriorPoint( S )
% Returns the most anterior point on a quadric surface
%
% Syntax:
%   X = quadric.mostAnteriorPoint( S )
%
% Description:
%   This is quite a kludge. I need to find the most anterior point on a
%   quadric surface (typically the cornea). This is hard to do in the
%   setting of a rotated quadric. So, I brute force it here.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   X                     - 1x3 vector containing the coordinates of the
%                           anterior-most point on the surface.
%

myFun = @(x) ellipseAreaAtX(S,x);

Scenter = quadric.center(S);
Sradii = quadric.radii(S);

p1 = fminbnd(myFun,Scenter(1), max(Sradii));

[~, c ] = myFun(p1);

X = [p1, c];

    
end


function [area, center] = ellipseAreaAtX(S,x)


% Obtain the variables for the quadric
[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVars(S);

% Obtain the elli        
a = B;
b = F;
c = C;
d = (D*x + H);
f = (E*x + I);
g = A*x.^2 + 2*G*x + K;


N = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g);
D = b^2-a*c;
S = realsqrt((a-c)^2+4*b^2);

check =  N/(D*(S-(a+c))) ;

if check < 0 
    area = nan;
    center = nan;
    return
end

transparent = ellipse_ex2transparent(ellipse_im2ex([a 2*b c 2*d 2*f g]));
area = transparent(3);
center = [transparent(1) transparent(2)];

end