function X = mostAnteriorPoint( S )
% Returns the most anterior point on an ellipsoidal surface
%
% Syntax:
%   X = quadric.mostAnteriorPoint( S )
%
% Description:
%   This is quite a kludge. I need to find the most anterior point on an
%   ellipsoidal surface (typically the cornea). This is hard to do in the
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


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Test here that the passed quadric is an ellipsoid
if ~strcmp(quadric.classify( S ),'ellipsoid')
    error('quadric:mostAnteriorPoint','Can only find the most anterior point for an ellipsoid');
end

% An anonymous function that returns the area of an ellipse that is in the
% cross-sectional plane at a given position along the x axis.
myFun = @(x) ellipseAreaAtX(S,x);

% Obtain the center and radii of the ellipsoid
Scenter = quadric.center(S);
Sradii = quadric.radii(S);

% Perform a search across the anonymous function, bound by the center and
% the largest radius. We find where the ellipse has an area that converges
% upon zero, and thus is the most anterior point.
p1 = fminbnd(myFun, Scenter(1), max(Sradii));

% Call the area function and now retain the [y z] coordinates of the
% ellipse center, which is the most anterior point
[~, c ] = myFun(p1);

% Assemble the 3D coordinate and return
X = [p1, c];

    
end


% Local function to return ellipse area
function [area, center] = ellipseAreaAtX(S,x)


% Obtain the variables for the quadric
[A, B, C, D, E, F, G, H, I, K] = quadric.matrixToVars(S);

% Map the ellipsoid parameters into the implicit parameters of an ellipse
% at the specfied x-axis position
a = B;
b = F;
c = C;
d = (D*x + H);
f = (E*x + I);
g = A*x.^2 + 2*G*x + K;

% Check that the ellipse is well formed
N = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g);
D = b^2-a*c;
S = realsqrt((a-c)^2+4*b^2);

check =  N/(D*(S-(a+c))) ;

if check < 0 
    area = nan;
    center = nan;
    return
end

% Convert the ellipse parameters from implicit to transparent format, and
% extract the area and ellipse center
transparent = ellipse_ex2transparent(ellipse_im2ex([a 2*b c 2*d 2*f g]));
area = transparent(3);
center = [transparent(1) transparent(2)];

end