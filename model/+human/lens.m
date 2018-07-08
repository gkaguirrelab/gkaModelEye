function lens = lens( eye )

% The front and back surfaces of the lens are modeled as
% hyperbolas, and as a set of ellipsoidal surfaces within the body
% of the lens with a gradient of refractive indices. All values
% taken from Atchison 2006.
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
lens.front.R = 11.48;
lens.front.Q = -5;
a = lens.front.R * sqrt(abs( 1 / (lens.front.Q - 1 ) )) * sign(lens.front.Q);
b = lens.front.R / (lens.front.Q - 1 );
radii(1) = abs(b);
radii(2:3) = abs(a);

S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[eye.pupil.center(1)+radii(1) 0 0]);
lens.front.S = quadric.matrixToVec(S);
lens.front.boundingBox = [-5.14 -3.7 -4 4 -4 4];
lens.front.side = 1;


lens.back.R = -5.9;
lens.back.Q = -2;
a = lens.back.R * sqrt(abs( 1 / (lens.back.Q - 1 ) )) * sign(lens.back.Q);
b = lens.back.R / (lens.back.Q - 1 );
radii(1) = abs(b);
radii(2:3) = abs(a);

S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[-7.3-radii(1) 0 0]);
lens.back.S = quadric.matrixToVec(S);
lens.back.boundingBox = [-7.3 -5.14 -4 4 -4 4];
lens.back.side = -1;

% I specify the location of a single nodal point to support
% calculation of the visual axis. The nodal point is placed at a
% depth of 7.2 mm, which is mid point of the nodal points specified
% in the Gullstrand-Emsley simplified schematic eye
lens.nodalPoint = [-7.2 0 0];

end

