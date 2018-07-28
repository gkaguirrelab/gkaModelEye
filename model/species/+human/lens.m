function lens = lens( eye, accomodation )

% Currently only supports odd number of shells.
nShells = 11;
startShell = 2;
endShell = 11;

% Initialize the components of the optical system
lens.S = [];
lens.boundingBox = [];
lens.side = [];
lens.mustIntersect = [];
lens.index = [];
lens.label = {};
lens.plot.color = {};

% Obtain the core and edge refractive indices
nEdge = returnRefractiveIndex( 'lens.edge', eye.meta.spectralDomain );
nCore = returnRefractiveIndex( 'lens.core', eye.meta.spectralDomain );

% This is the position (on the optical axis) of the point in the lens with
% the maximal refractive index.
lensCenter = -5.4;

% The thickness of the back and front of the lens.
lensThickBack = 2.16;
lensThickFront = 1.44;


%% Back lens surface
% The back surface of the lens is modeled as a hyperbola. Values taken
% from Atchison 2006.
% To convert R and Q to radii of a hyperbola:
%   R = b^2/a
%	Q = (a^2 / b^2) + 1
% Therefore, given R and Q, we can obtain a and b, which correspond
% to the radii of the ellipsoid model, with a corresponding to the
% axial dimension, and b to the horizontal and verical dimensions.
% Checking my algebra here:
%{
            syms a b R Q
            eqn1 = R == a^2/b;
            eqn2 = Q == (a^2 / b^2) + 1;
            solution = solve([eqn1, eqn2]);
            solution.a
            solution.b
%}
R = -5.9; Q = -2;
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
radii(1) = abs(b); radii(2:3) = abs(a);

% Build the quadric
S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[lensCenter-lensThickBack-radii(1) 0 0]);
boundingBox = [lensCenter-lensThickBack lensCenter -5 5 -5 5];

% Add to the optical system structure
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; -1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.index = [lens.index; nEdge];
lens.label = [lens.label; {'lens.back'}];
lens.plot.color = [lens.plot.color; {'red'}];

%% Back gradient shells
boundingBox = [lensCenter-lensThickBack lensCenter -5 5 -5 5];
nLensVals = linspace(nEdge,nCore,nShells*2+1);
for ii = startShell:endShell
    nBackQuadric = [-0.010073731138546 0 0 0.0; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0 0 0 1.418000000000000];
    S = nBackQuadric;
    S(end,end)=S(end,end)-nLensVals(ii*2);
    S = quadric.translate(S,[lensCenter 0 0]);
    
    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; 1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii*2+1)];
    lens.label = [lens.label; {sprintf('lens.back.shell_n=%0.3f',nLensVals(ii*2))}];
    lens.plot.color = [lens.plot.color; [0.5 (nLensVals(ii*2)-nEdge)/(nCore-nEdge)/2+0.5 0.5]];
end
% Force the index of the center shell to be equal to the lens core
% refractive index. The value can be something other than this when
% endShell shell is not equal to the nShells.
lens.index(end) = nCore;


%% Front gradient shells
boundingBox = [lensCenter lensCenter+lensThickFront -5 5 -5 5];
nLensVals = linspace(nEdge,nCore,nShells*2+1);
for ii = endShell:-1:startShell
    nFrontQuadric = [0.022665895061728 0 0 0; 0 0.002039930555556 0 0; 0 0 0.002039930555556 0; 0 0 0 1.371000000000000];
    S = nFrontQuadric;
    S(end,end)=nLensVals(ii*2)-nCore;
    S = quadric.translate(S,[lensCenter-0.065277777777778 0 0]);
    
    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; 1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii*2-1)];
    lens.label = [lens.label; {sprintf('lens.front.shell_n=%0.3f',nLensVals(ii*2))}];
    lens.plot.color = [lens.plot.color; [0.5 (nLensVals(ii*2)-nEdge)/(nCore-nEdge)/2+0.5 0.5]];
end
% Force the index of the space between the last front shell and the front
% lens surface to be equal to the lens edge refractive index. The value can
% be something other than this when the start shell is something other than
% one.
lens.index(end) = nEdge;


%% Front lens surface
R = 11.48;
Q = -5;
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
radii(1) = abs(b);
radii(2:3) = abs(a);

% Build the quadric
S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[eye.pupil.center(1)+radii(1) 0 0]);
c = quadric.center(S); r = quadric.radii(S);
boundingBox = [lensCenter c(1)-r(1) -5 5 -5 5];

% Add to the optical system structure. No refractive index added with this
% surface, as this is the last surface of this set.
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; 1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.label = [lens.label; {'lens.front'}];
lens.plot.color = [lens.plot.color; {'red'}];


end

