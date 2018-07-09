function lens = lens( eye )

nShells = 5;

% Initialize the components of the optical system
lens.S = [];
lens.boundingBox = [];
lens.side = [];
lens.mustIntersect = [];
lens.index = [];
lens.label = {};

% Obtain the core and edge refractive indices
nEdge = returnRefractiveIndex( 'lens.edge', eye.meta.spectralDomain );
nCore = returnRefractiveIndex( 'lens.core', eye.meta.spectralDomain );

% This is the position (on the optical axis) of the point in the lens with
% the maximal refractive index.
lensCenter = -5.4;
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
boundingBox = [lensCenter-lensThickBack lensCenter -4 4 -4 4];

% Add to the optical system structure
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; -1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.index = [lens.index; nEdge];
lens.label = {lens.label{:}; 'lens.back'};


%% Back gradient shells
boundingBox = [lensCenter-lensThickBack lensCenter -4 4 -4 4];
nLensVals = linspace(nEdge,nCore,nShells+1);
for ii = 1:nShells
    nBackQuadric = [-0.010073731138546 0 0 0.0; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0 0 0 1.418000000000000];
    S = nBackQuadric;
    S(end,end)=S(end,end)-nLensVals(ii)+1;
    S=S./S(end,end);
    S = quadric.translate(S,[lensCenter 0 0]);
    
    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; 1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii+1)];
    lens.label = [lens.label; {sprintf('lens.back.shell%02d',ii)}];
end


%% Front gradient shells
boundingBox = [lensCenter lensCenter+lensThickFront -4 4 -4 4];
nLensVals = linspace(nCore,nEdge,nShells+11);
for ii = 1:nShells
    nFrontQuadric = [0.022665895061728 0 0 0; 0 0.002039930555556 0 0; 0 0 0.002039930555556 0; 0 0 0 1.371000000000000];
    S = nFrontQuadric;
    S(end,end)=nLensVals(ii)-nCore;
    S = quadric.translate(S,[lensCenter-0.065277777777778 0 0]);
    S=S./S(end,end);

    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; 1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii+1)];
    lens.label = [lens.label; {sprintf('lens.front.shell%02d',ii)}];
end


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
boundingBox = [lensCenter lensCenter+lensThickFront -4 4 -4 4];

% Add to the optical system structure
% No refractive index added with this surface, as this is the last surface
% of this set.
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; 1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.label = [lens.label; {'lens.front'}];




% I specify the location of a single nodal point to support
% calculation of the visual axis. The nodal point is placed at a
% depth of 7.2 mm, which is mid point of the nodal points specified
% in the Gullstrand-Emsley simplified schematic eye
lens.nodalPoint = [-7.2 0 0];

end

