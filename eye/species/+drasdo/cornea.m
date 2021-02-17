function cornea = cornea( eye )
% Returns the cornea sub-field of an eye model structure
%
% Syntax:
%  cornea = human.cornea( eye )
%
% Description:
%   
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   cornea                - Structure.
%





%% Front corneal surface
% The Drasdo & Fowler model describes a corneal front surface that has a
% radius of curvature of 7.8 mm, and an eccentricity of 0.5. This code
% obtains the semi-radi for this ellipse:
%{
    eqn1 = e == sqrt(1 - (b^2/a^2));
    eqn2 = r == (b^2)/a;
    eqn3 = r == 7.8;
    eqn4 = e == 0.5;
    sol=solve([eqn1, eqn2, eqn3, eqn4]);
    radii = [eval(sol.a(2)), eval(sol.b(2)), eval(sol.b(2))];
%}

% Corneal radii
radii = [10.4000    9.0067    9.0067];

% Create the quadric
S = quadric.scale(quadric.unitSphere,radii);

% We set the center of the cornea front surface ellipsoid so that the axial
% apex (prior to rotation) is at position [0, 0, 0]
S = quadric.translate(S,[-radii(1) 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store these values
cornea.front.S = quadric.matrixToVec(S);
cornea.front.side = 1;
cornea.front.boundingBox=[-4 X(1) -9 9 -9 9];
cornea.front.center=[-radii(1) 0 0];

% Create a back corneal surface, which is just the front translated a tiny
% bit posteriorly
delta = -0.01;
S = quadric.translate(S,[delta 0 0]);
cornea.back.S = quadric.matrixToVec(S);
cornea.back.side = 1;
cornea.back.boundingBox=[-4-delta X(1)-delta -9 9 -9 9];
cornea.back.center=[-radii(1)-delta 0 0];



%% Store the combined corneal surfaces
cornea.S = [cornea.back.S; cornea.front.S];
cornea.boundingBox = [cornea.back.boundingBox; cornea.front.boundingBox];
cornea.side = [1; 1];
cornea.mustIntersect = [1; 1];
cornea.index = 1.336;
cornea.label = [{'cornea.back'};{'cornea.front'}];
cornea.plot.color = [{[0.5 0.5 0.75]};{[0.5 0.5 0.75]}];


end