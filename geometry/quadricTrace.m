clear opticalSystem
clear intersectionCoords
clear rayPath

opticalSystem(1,:)=[nan(1,10) nan nan(1,6) nan 1.357];


% %% Back lens surface
S = quadric.unitTwoSheetHyperboloid();
S = quadric.scale(S,[ 1.9667 3.4064 3.4064]);
S = quadric.translate(S,[-7.3-1.9667 0 0]);
boundingBox=[-7.4 -5.14 -4 4 -4 4];
opticalSystem(end+1,:)=[quadric.matrixToVec(S) -1 boundingBox 1 1.371];

%% Back lens gradient shells
nLensVals = linspace(1.371,1.418,7);
for ii = 2:length(nLensVals)-1
    nBackQuadric = [-0.010073731138546 0 0 0.0; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0 0 0 1.418000000000000];
    S = nBackQuadric;
    S(end,end)=S(end,end)-nLensVals(ii);
    S=S./S(end,end);
    S = quadric.translate(S,[-5.14 0 0]);
    boundingBox=[-7.3 -5.14 -4 4 -4 4];
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 0 nLensVals(ii+1)];
end

nLensVals = linspace(1.418,1.371,7);
for ii = 2:length(nLensVals)-1
    nFrontQuadric = [0.022665895061728 0 0 0; 0 0.002039930555556 0 0; 0 0 0.002039930555556 0; 0 0 0 1.371000000000000];
    S = nFrontQuadric;
    S(end,end)=nLensVals(ii)-1.418;
    S = quadric.translate(S,[-5.14-0.065277777777778 0 0]);
    S=S./S(end,end);
    boundingBox=[-5.14 -3.5 -4 4 -4 4];
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 0 nLensVals(ii)];
end

%% Front lens surface
S = quadric.unitTwoSheetHyperboloid();
S = quadric.scale(S,[ 1.9133 4.6867 4.6867]);
S = quadric.translate(S,[-3.7+1.9133 0 0]);
boundingBox=[-5.14 -3.5 -4 4 -4 4];
opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1 1.3335];

S = quadric.unitSphere();
S = quadric.scale(S,[ 13.7716  9.3027  9.3027]);
S = quadric.translate(S,[-14.3216 0 0]);
boundingBox=[-4 0 -8 8 -8 8];
opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1 1.3747];

S = quadric.unitSphere();
S = quadric.scale(S,[ 14.2600  10.4300 10.2700    ]);
S = quadric.translate(S,[-14.2600 0 0]);
boundingBox=[-4 0 -8 8 -8 8];
opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1 1.0];

% Plot the optical system
plotOpticalSystem('opticalSystem',opticalSystem,'addLighting',true);

% Trace a bundle of rays starting from a retinal point

for angle = -10:2.5:10

p = [-23.58;0;0];
u = [1;tand(angle);0];
u = u./sqrt(sum(u.^2));
R = [p, u];

[outputRay, rayPath] = rayTraceQuadrics(R, opticalSystem);
plotOpticalSystem('newFigure',false,'outputRay',outputRay,'rayPath',rayPath);
end

