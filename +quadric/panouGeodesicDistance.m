function [distance,startAngle,endAngle,geodeticPathCoords] = panouGeodesicDistance(S,X0,X1)
% Find the geodesic distance between two points on a tri-axial ellipsoid
%
% Syntax:
%  [distance,startAngle,endAngle,geodeticPathCoords] = panouGeodesicDistance(S,X0,X1)
%
% Description:
%   Returns the geodesic distance between two points on the tri-axial
%   ellipsoidal surface. The code was provided by Georgios Panou, and is
%   based upon the approach described in:
%
%       Panou, G. "The geodesic boundary value problem and its solution on
%       a triaxial ellipsoid." Journal of Geodetic Science 3.3 (2013):
%       240-249.
%
% Inputs:
%   S                     - 4x4 quadric surface matrix
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%
% Outputs:
%   distance              - Scalar. Distance of the geodetic between the
%                           two points.
%   startAngle, endAngle  - 
%
% Examples:
%{
    % Numeric example provided by G. Panou in the original code
    S = quadric.scale(quadric.unitSphere,[0.015, 0.010, 0.009]);
    X0 = quadric.ellipsoidalGeoToCart( [5 5 0], S );
    X1 = quadric.ellipsoidalGeoToCart( [60 120 0], S );
    [distance,startAngle,endAngle,geodeticPathCoords] = quadric.panouGeodesicDistance(S,X0,X1);
    assert( max(abs(s - 0.0259)) < 1e-3 );

        F = quadric.vecToFunc(quadric.matrixToVec(S));
    figure
    [xx, yy, zz] = meshgrid( linspace(-2,2,100), linspace(-2,2,100), linspace(-2,2,100));
    vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);
    p = patch(vertices);
    p.FaceColor =[.9 .9 .9];
    p.EdgeColor = 'none';
    alpha(.8);
    daspect([1 1 1])
    view([40 30]); 
    axis tight
    axis equal
    camlight
    lighting gouraud
    hold on
    % Plot lines of varying beta
    for beta = -90:10:90
        coords =[];
        for omega = -180:3:180
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        plot3(coords(:,1),coords(:,2),coords(:,3),'-b');
    end
    fprintf('Lines of constant beta in blue\n');
    % Plot lines of varying omega
    for omega = -180:10:180
        coords =[];
        for beta = -90:3:90
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        plot3(coords(:,1),coords(:,2),coords(:,3),'-g');
    end
    fprintf('Lines of constant omega in green\n');
    plot3(geodeticPathCoords(:,1),geodeticPathCoords(:,2),geodeticPathCoords(:,3),'-r');
%}


% number of subspaces
nSubspaces = 20000;

% desired accuracy (difference of beta1 in radians). For very small bodies,
% may be larger, e.g. 10^(-9)
epsilonAccuracy = 10^(-12);

% resolution of the returned path coords
pathResolution = 20;

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Obtain the radii of the quadric surface and distribute the values. We
% adopt the canonical order of a => b => c, but the order returned after
% alignment of the axes is a <= b <= c. This is why c is mapped to the
% first value in the radii.
radii = quadric.radii(quadric.alignAxes(S));
a=radii(3);b=radii(2);c=radii(1);

% Obtain the ellipsoidal geodetic coordinates for the two points, and
% convert to radians
G0 = deg2rad(quadric.cartToEllipsoidalGeo( X0, S ));
G1 = deg2rad(quadric.cartToEllipsoidalGeo( X1, S ));

if G0(2) == G1(2)
    [nIterations,LiouvilleConstant,startAngle,endAngle,distance,geodeticPathCoords] = ...
        Geodesics_dbeta(a,b,c,G0(1),G0(2),G1(1),G1(2),nSubspaces,epsilonAccuracy);
else
    [nIterations,LiouvilleConstant,startAngle,endAngle,distance,geodeticPathCoords] = ...
        Geodesics(a,b,c,G0(1),G0(2),G1(1),G1(2),nSubspaces,epsilonAccuracy);
end

% Prepare a subset of the geodeticPathCoords to return
sampleBase = round(linspace(1,size(geodeticPathCoords,1),pathResolution));
geodeticPathCoords=geodeticPathCoords(sampleBase',:);

% Re-arrange the order of the dimensions of the geodeticPathCoords so that
% they reflect the axis order of the passed quadric surface
geodeticPathCoords = geodeticPathCoords(:,[3 2 1]);

% Rotate and translate the geodetic path coords so that they correspond to
% the surface of the original quadric
% Now counter-rotate and then counter-translate the X coordinate
[~, rotMat] = quadric.alignAxes(S);
for ii = 1:pathResolution
    geodeticPathCoords(ii,:) = (geodeticPathCoords(ii,:)*rotMat');
end

origCenter = quadric.center(S);
geodeticPathCoords = geodeticPathCoords+origCenter';

end % panouGeodesicDistance




function [n,c,a0,a1,s,geodeticPathCoords] = Geodesics_dbeta(ax,ay,b,beta0,lambda0,beta1,lambda1,Sub,epsilon)

%----- Ellipsoid -----
hx = sqrt(ax^2-b^2);
hy = sqrt(ay^2-b^2);
he = sqrt(ax^2-ay^2);

[beta,lambda] = meshgrid(linspace(0,2*pi,25),linspace(0,2*pi,25));
x = ax.*sqrt(cos(beta).^2+((he/hx)^2).*sin(beta).^2).*cos(lambda);
y = ay.*cos(beta).*sin(lambda);
z = b.*sin(beta).*sqrt(1-((he/hx)^2).*cos(lambda).^2);
mesh(x,y,z)
hold on
%----------

omega = acos(sin(lambda0)*sin(lambda1)+cos(lambda0)*cos(lambda1)*cos(beta1-beta0));
Alpha = acos((sin(lambda1)-sin(lambda0)*cos(omega))/(cos(lambda0)*sin(omega)));

dlambda0 = cos(lambda0)*cot(Alpha);

ddlambda0 = 0;

n=1;
while n<100
    n;
    
    Dlambda0 = dlambda0;
    Ddlambda0 = ddlambda0;
    
    R = RK4_dbeta(beta0, beta1, Sub, [lambda0,Dlambda0+Ddlambda0,0,1], ax, ay, b);
    
    beta = R(:,1);
    lambda = R(:,2);
    dlambda = R(:,3);
    plambda = R(:,4);
    pdlambda = R(:,5);
    
    lambda10 = lambda(end);
    plambda10 = plambda(end);
    
    dlambda0 = dlambda(1);
    ddlambda0 = (lambda1-lambda10)/plambda10;
    
    %Check
    n = n+1;
    if abs(lambda1-lambda10) < epsilon
        break
    end
end
%----------

x = ax.*sqrt(cos(beta).^2+((he/hx)^2).*sin(beta).^2).*cos(lambda);
y = ay.*cos(beta).*sin(lambda);
z = b.*sin(beta).*sqrt(1-((he/hx)^2).*cos(lambda).^2);

geodeticPathCoords = [x y z];

plot3(x,y,z)

B = (ay^2.*sin(beta).^2+b^2.*cos(beta).^2)./(hx^2-hy^2.*sin(beta).^2);
L = (ax^2.*sin(lambda).^2+ay^2.*cos(lambda).^2)./(hx^2-he^2.*cos(lambda).^2);

E = B.*(hy^2.*cos(beta).^2+he^2.*sin(lambda).^2);
G = L.*(hy^2.*cos(beta).^2+he^2.*sin(lambda).^2);

%----- Liouville's constant, Azimuths -----

alpha = acot(sqrt(L./B).*(dlambda));

for i = 1:1:length(alpha)
    if alpha(i) < 0
        alpha(i) = alpha(i) + pi;
    end
end

c2 = (cos(beta).^2+((he/hx)^2).*sin(beta).^2).*(cos(alpha).^2)+((he/hx)^2).*(cos(lambda).^2).*(sin(alpha).^2);
c = sqrt(c2);

c = c(Sub/2);

alpha = alpha.*(180/pi);

a0 = alpha(1);

a1 = alpha(end);

%----------

%----- Geodesic distance -----

h = (beta1-beta0)/Sub;

Fv = sqrt(E+G.*(dlambda).^2);

S = [1 zeros(1,Sub-1) 1];
S(2:2:Sub) = 4*ones(1,Sub/2);
S(3:2:Sub-1) = 2*ones(1,Sub/2-1);

s = (h/3)*sum(S*Fv);

end


function [n,c,a0,a1,s,geodeticPathCoords] = Geodesics(ax,ay,b,beta0,lambda0,beta1,lambda1,Sub,epsilon)

%----- Ellipsoid -----
hx = sqrt(ax^2-b^2);
hy = sqrt(ay^2-b^2);
he = sqrt(ax^2-ay^2);

[beta,lambda] = meshgrid(linspace(0,2*pi,25),linspace(0,2*pi,25));
x = ax.*sqrt(cos(beta).^2+((he/hx)^2).*sin(beta).^2).*cos(lambda);
y = ay.*cos(beta).*sin(lambda);
z = b.*sin(beta).*sqrt(1-((he/hx)^2).*cos(lambda).^2);
mesh(x,y,z)
hold on
%----------

omega = acos(sin(beta0)*sin(beta1)+cos(beta0)*cos(beta1)*cos(lambda1-lambda0));
Alpha = acos((sin(beta1)-sin(beta0)*cos(omega))/(cos(beta0)*sin(omega)));

dbeta0 = cos(beta0)*cot(Alpha);

ddbeta0 = 0;

n=1;
while n<100
    n;
    
    Dbeta0 = dbeta0;
    Ddbeta0 = ddbeta0;
    
    R = RK4(lambda0, lambda1, Sub, [beta0,Dbeta0+Ddbeta0,0,1], ax, ay, b);
    
    lambda = R(:,1);
    beta = R(:,2);
    dbeta = R(:,3);
    pbeta = R(:,4);
    pdbeta = R(:,5);
    
    beta10 = beta(end);
    pbeta10 = pbeta(end);
    
    dbeta0 = dbeta(1);
    ddbeta0 = (beta1-beta10)/pbeta10;
    
    %Check
    n = n+1;
    if abs(beta1-beta10) < epsilon
        break
    end
end
%----------

x = ax.*sqrt(cos(beta).^2+((he/hx)^2).*sin(beta).^2).*cos(lambda);
y = ay.*cos(beta).*sin(lambda);
z = b.*sin(beta).*sqrt(1-((he/hx)^2).*cos(lambda).^2);
plot3(x,y,z)
geodeticPathCoords = [x y z];


B = (ay^2.*sin(beta).^2+b^2.*cos(beta).^2)./(hx^2-hy^2.*sin(beta).^2);
L = (ax^2.*sin(lambda).^2+ay^2.*cos(lambda).^2)./(hx^2-he^2.*cos(lambda).^2);

E = B.*(hy^2.*cos(beta).^2+he^2.*sin(lambda).^2);
G = L.*(hy^2.*cos(beta).^2+he^2.*sin(lambda).^2);

%----- Liouville's constant, Azimuths -----

alpha = acot(sqrt(B./L).*(dbeta));

for i = 1:1:length(alpha)
    if alpha(i) < 0
        alpha(i) = alpha(i) + pi;
    end
end

c2 = (cos(beta).^2+((he/hx)^2).*sin(beta).^2).*(sin(alpha).^2)+((he/hx)^2).*(cos(lambda).^2).*(cos(alpha).^2);
c = sqrt(c2);

c = c(Sub/2);

alpha = alpha.*(180/pi);

a0 = alpha(1);

a1 = alpha(end);

%----------

%----- Geodesic distance -----

h = (lambda1-lambda0)/Sub;

Fv = sqrt(E.*(dbeta).^2+G);

S = [1 zeros(1,Sub-1) 1];
S(2:2:Sub) = 4*ones(1,Sub/2);
S(3:2:Sub-1) = 2*ones(1,Sub/2-1);

s = (h/3)*sum(S*Fv);
end


function dx = F_dbeta(beta,x,ax,ay,b)
dx = zeros(4,1);

%----- Ellipsoid -----
hx = sqrt(ax^2-b^2);
hy = sqrt(ay^2-b^2);
he = sqrt(ax^2-ay^2);
%----------

B = (ay^2.*sin(beta).^2+b^2.*cos(beta).^2)./(hx^2-hy^2.*sin(beta).^2);
L = (ax^2.*sin(x(1)).^2+ay^2.*cos(x(1)).^2)./(hx^2-he^2.*cos(x(1)).^2);

Bp = (ax^2.*hy^2.*sin(2.*beta))./((hx^2-hy^2.*sin(beta).^2).^2);
Lp = -(b^2.*he^2.*sin(2.*x(1)))./((hx^2-he^2.*cos(x(1)).^2).^2);

Bpp = (2.*ax^2.*hy^4.*sin(2.*beta).^2)./((hx^2-hy^2.*sin(beta).^2).^3)+(2.*ax^2.*hy^2.*cos(2.*beta))./((hx^2-hy^2.*sin(beta).^2).^2);
Lpp = (2.*b^2.*he^4.*sin(2.*x(1)).^2)./((hx^2-he^2.*cos(x(1)).^2).^3)-(2.*b^2.*he^2.*sin(2.*x(1)))./((hx^2-he^2.*cos(x(1)).^2).^2);

%----------
E = B.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2);

Eb = Bp.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2)-B.*hy^2.*sin(2.*beta);
El = B.*he^2.*sin(2.*x(1));

Ebb = Bpp.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2)-2.*Bp.*hy^2.*sin(2.*beta)-2.*B.*hy^2.*cos(2.*beta);
Ell = 2.*B.*he^2.*cos(2.*x(1));

Ebl = Bp.*he^2.*sin(2.*x(1));
Elb = Bp.*he^2.*sin(2.*x(1));

%----------
G = L.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2);

Gl = Lp.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2)+L.*he^2.*sin(2.*x(1));
Gb = -L.*hy^2.*sin(2.*beta);

Gll = Lpp.*(hy^2.*cos(beta).^2+he^2.*sin(x(1)).^2)+2.*Lp.*he^2.*sin(2.*x(1))+2.*L.*he^2.*cos(2.*x(1));
Gbb = -2.*L.*hy^2.*cos(2.*beta);

Glb = -Lp.*hy^2.*sin(2.*beta);
Gbl = -Lp.*hy^2.*sin(2.*beta);

%----------

dx(1) = x(2);


q3 = -(1/2)*(Gb/E);

q2 = (El/E)-(1/2)*(Gl/G);

q1 = (1/2)*(Eb/E)-(Gb/G);

q0 = (1/2)*(El/G);

dx(2) = q3*x(2)^3+q2*x(2)^2+q1*x(2)+q0;


dx(3) = x(4);


q33 = -(1/2)*((E*Gbl-El*Gb)/(E^2));

q22 = ((E*Ell-El*El)/(E^2))-(1/2)*((G*Gll-Gl*Gl)/(G^2));

q11 = (1/2)*((E*Ebl-Eb*El)/(E^2))-((G*Gbl-Gb*Gl)/(G^2));

q00 = (1/2)*((G*Ell-El*Gl)/(G^2));

dx(4) = (q33*x(2)^3+q22*x(2)^2+q11*x(2)+q00)*x(3)+(3*q3*x(2)^2+2*q2*x(2)+q1)*x(4);

end



function dx = F(lambda,x,ax,ay,b)
dx = zeros(4,1);

%----- Ellipsoid -----
hx = sqrt(ax^2-b^2);
hy = sqrt(ay^2-b^2);
he = sqrt(ax^2-ay^2);
%----------

B = (ay^2.*sin(x(1)).^2+b^2.*cos(x(1)).^2)./(hx^2-hy^2.*sin(x(1)).^2);
L = (ax^2.*sin(lambda).^2+ay^2.*cos(lambda).^2)./(hx^2-he^2.*cos(lambda).^2);

Bp = (ax^2.*hy^2.*sin(2.*x(1)))./((hx^2-hy^2.*sin(x(1)).^2).^2);
Lp = -(b^2.*he^2.*sin(2.*lambda))./((hx^2-he^2.*cos(lambda).^2).^2);

Bpp = (2.*ax^2.*hy^4.*sin(2.*x(1)).^2)./((hx^2-hy^2.*sin(x(1)).^2).^3)+(2.*ax^2.*hy^2.*cos(2.*x(1)))./((hx^2-hy^2.*sin(x(1)).^2).^2);
Lpp = (2.*b^2.*he^4.*sin(2.*lambda).^2)./((hx^2-he^2.*cos(lambda).^2).^3)-(2.*b^2.*he^2.*sin(2.*lambda))./((hx^2-he^2.*cos(lambda).^2).^2);

%----------
E = B.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2);

Eb = Bp.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2)-B.*hy^2.*sin(2.*x(1));
El = B.*he^2.*sin(2.*lambda);

Ebb = Bpp.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2)-2.*Bp.*hy^2.*sin(2.*x(1))-2.*B.*hy^2.*cos(2.*x(1));
Ell = 2.*B.*he^2.*cos(2.*lambda);

Ebl = Bp.*he^2.*sin(2.*lambda);
Elb = Bp.*he^2.*sin(2.*lambda);

%----------
G = L.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2);

Gl = Lp.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2)+L.*he^2.*sin(2.*lambda);
Gb = -L.*hy^2.*sin(2.*x(1));

Gll = Lpp.*(hy^2.*cos(x(1)).^2+he^2.*sin(lambda).^2)+2.*Lp.*he^2.*sin(2.*lambda)+2.*L.*he^2.*cos(2.*lambda);
Gbb = -2.*L.*hy^2.*cos(2.*x(1));

Glb = -Lp.*hy^2.*sin(2.*x(1));
Gbl = -Lp.*hy^2.*sin(2.*x(1));

%----------

dx(1) = x(2);


p3 = -(1/2)*(El/G);

p2 = (Gb/G)-(1/2)*(Eb/E);

p1 = (1/2)*(Gl/G)-(El/E);

p0 = (1/2)*(Gb/E);

dx(2) = p3*x(2)^3+p2*x(2)^2+p1*x(2)+p0;


dx(3) = x(4);


p33 = -(1/2)*((G*Elb-El*Gb)/(G^2));

p22 = ((G*Gbb-Gb*Gb)/(G^2))-(1/2)*((E*Ebb-Eb*Eb)/(E^2));

p11 = (1/2)*((G*Glb-Gb*Gl)/(G^2))-((E*Elb-Eb*El)/(E^2));

p00 = (1/2)*((E*Gbb-Eb*Gb)/(E^2));

dx(4) = (p33*x(2)^3+p22*x(2)^2+p11*x(2)+p00)*x(3)+(3*p3*x(2)^2+2*p2*x(2)+p1)*x(4);
end


function R = RK4_dbeta(beta0, beta1, Sub, x0, ax, ay, b)
%----- Runge–Kutta_4 method -----

h = (beta1-beta0)/Sub;      %step size
beta(1) = beta0;
w(:,1) = x0;                    %initial conditions

for i = 1:Sub
    k1 = h*F_dbeta(beta(i), w(:,i), ax,ay,b);
    k2 = h*F_dbeta(beta(i)+0.5*h, w(:,i)+0.5*k1, ax,ay,b);
    k3 = h*F_dbeta(beta(i)+0.5*h, w(:,i)+0.5*k2, ax,ay,b);
    k4 = h*F_dbeta(beta(i)+h, w(:,i)+k3, ax,ay,b);
    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    beta(i+1) = beta0 + i*h;
end

R = [beta' w'];

end


function R = RK4(lambda0, lambda1, Sub, x0, ax, ay, b)
%----- Runge–Kutta_4 method -----

h = (lambda1-lambda0)/Sub;      %step size
lambda(1) = lambda0;
w(:,1) = x0;                    %initial conditions

for i = 1:Sub
    k1 = h*F(lambda(i), w(:,i), ax,ay,b);
    k2 = h*F(lambda(i)+0.5*h, w(:,i)+0.5*k1, ax,ay,b);
    k3 = h*F(lambda(i)+0.5*h, w(:,i)+0.5*k2, ax,ay,b);
    k4 = h*F(lambda(i)+h, w(:,i)+k3, ax,ay,b);
    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    lambda(i+1) = lambda0 + i*h;
end

R = [lambda' w'];
end