clear opticalSystem
clear intersectionCoords


opticalSystem(:,1)=[nan(1,10) nan nan(1,6) 1.3370];

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
boundingBox=[-8 -3.5 -4 4 -4 4];
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 1 boundingBox 1.4];

S = quadric.unitTwoSheetHyperboloid();
S = quadric.scale(S,[ 1.9133 4.6867 4.6867]);
S = quadric.translate(S,[-2.0117 0 0]);
boundingBox=[-8 -3.5 -4 4 -4 4];
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 boundingBox 1.3];

S = quadric.unitSphere();
S = quadric.scale(S,[ 13.7716  9.3027  9.3027]);
S = quadric.translate(S,[-14.3216 0 0]);
boundingBox=[-4 0 -8 8 -8 8];
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 boundingBox 1.3747];

S = quadric.unitSphere();
S = quadric.scale(S,[ 14.2600  10.4300 10.2700    ]);
S = quadric.translate(S,[-14.2600 0 0]);
boundingBox=[-4 0 -8 8 -8 8];
opticalSystem(:,end+1)=[quadric.matrixToVec(S) 2 boundingBox 1.0];


p = [-3.925;0;0];
u = [1;tand(-15);tand(10)];
u = u./sqrt(sum(u.^2));
R = [p, u];
atan2(R(2,2),R(1,2));

figure

surfColors = {'',[0.5 0.5 0.5],[0.5 0.5 0.5],'blue','blue'};

rayPath(1,:) = p;

for ii=2:size(opticalSystem,2)
    % Extract components from optical system vector
    S = quadric.vecToMatrix(opticalSystem(1:10,ii));
    side = opticalSystem(11,ii);
    boundingBox = opticalSystem(12:17,ii);
    nRel = opticalSystem(18,ii-1)/opticalSystem(18,ii);

    % Compute the intersection, surface normal, and refracted rat
    X = quadric.intersectRay(S,R,side,boundingBox);
    N = quadric.surfaceNormal(S,X,side);
    R = quadric.refractRay(R,N,nRel);

    % Store the ray path
    rayPath(end+1,:) = X;
    
    % Obtain a function handle for the polynomial
    F = quadric.vecToFunc(quadric.matrixToVec(S));

    % Plot the surface
    plotSurface(F,boundingBox,surfColors{ii})
    hold on
        
end

% Plot the ray
rayPath(end+1,:) = R(:,1)+R(:,2)*5;

plot3(rayPath(:,1),rayPath(:,2),rayPath(:,3),'-r');

function plotSurface(F,boundingBox,surfColor)

[xx, yy, zz]=meshgrid( linspace(boundingBox(1),boundingBox(2),100),...
    linspace(boundingBox(3),boundingBox(4),100),...
    linspace(boundingBox(5),boundingBox(6),100));
    vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);

p = patch(vertices);
p.FaceColor = surfColor;
p.EdgeColor = 'none';
alpha(0.1);
daspect([1 1 1])
view(3); 
axis tight
axis equal
camlight 
lighting gouraud

end
