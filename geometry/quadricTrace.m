clear opticalSystem
opticalSystem(:,1)=[nan(1,10) nan 1.3370];

% Replicate Elagha
% S = quadric.unitSphere();
% S = quadric.scale(S,10);
% S = quadric.translate(S,[22; 0; 0]);
% opticalSystem(:,end+1)=[quadric.matrixToVec(S) 1 1.2];
% 
% S = quadric.unitSphere();
% S = quadric.scale(S,8);
% S = quadric.translate(S,[9; 0; 0]);
% opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 1];
% 
% S = quadric.unitSphere();
% S = quadric.scale(S,12);
% S = quadric.translate(S,[34; 0; 0]);
% opticalSystem(:,end+1)=[quadric.matrixToVec(S) 1 1.5];
% 
% S = quadric.unitSphere();
% S = quadric.scale(S,10);
% S = quadric.translate(S,[20; 0; 0]);
% opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 1.0];

S = quadric.unitTwoSheetHyperboloid();
S = quadric.scale(S,[ 1.9667 3.4064 3.4064]);
S = quadric.translate(S,[-9.4917 0 0]);
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 1 1.3747];

S = quadric.unitTwoSheetHyperboloid();
S = quadric.scale(S,[ 1.9133 4.6867 4.6867]);
S = quadric.translate(S,[-2.0117 0 0]);
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 1.3747];

S = quadric.unitSphere();
S = quadric.scale(S,[ 13.7716  9.3027  9.3027]);
S = quadric.translate(S,[-14.3216 0 0]);
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 1.3747];

S = quadric.unitSphere();
S = quadric.scale(S,[ 14.2600  10.4300 10.2700    ]);
S = quadric.translate(S,[-14.2600 0 0]);
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 1.0];


p = [-3.925;2;0];
u = [1;tand(-15);tand(10)];
u = u./sqrt(sum(u.^2));
R = [p, u];
atan2(R(2,2),R(1,2));

figure
    gv = linspace(15,-8,100); % adjust for appropriate domain
[xx, yy, zz]=meshgrid(gv, gv, gv);

surfColors = {'','blue','cyan','red','green'};

for ii=2:5
    S = quadric.vecToMatrix(opticalSystem(1:10,ii));
    side = opticalSystem(11,ii);
    X = quadric.intersectRay(S,R,side);
    N = quadric.surfaceNormal(S,X,side);
    nRel = opticalSystem(12,ii-1)/opticalSystem(12,ii);
    R = quadric.refractRay(R,N,nRel);

    % Obtain a function handle for the polynomial
    F = quadric.vecToFunc(quadric.matrixToVec(S));

    % Plot the surface
    vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);
    plotSurface(vertices,surfColors{ii})
    hold on
    axis equal

end


function plotSurface(vertices,surfColor)

p = patch(vertices);
%isonormals(x,y,z,v,p)
p.FaceColor = surfColor;
p.EdgeColor = 'none';
alpha(0.5);
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

end
