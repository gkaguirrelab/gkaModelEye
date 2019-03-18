function lens = lens( eye )
% Returns the lens sub-field of an eye model structure
%
% Syntax:
%  lens = human.lens( eye )
%
% Description:
%   A model of the crystalline lens is generated, expressed as a set of
%   quadric surfaces. The anterior and posterior surfaces of the lens are
%   modeled as one-half of a two-sheeted hyperboloid. The radii of the
%   surfaces, and their dependence upon age and the accommodative state of
%   the eye, are taken from:
%
%       Navarro, Rafael. "Adaptive model of the aging emmetropic eye and
%       its changes with accommodation." Journal of vision 14.13 (2014):
%       21-21.
%
%   The total lens thickness, and thickness of the anterior and posterior
%   portions of the lens, are also taken from Navarro 2014.
%
%   The interior of the lens is modeled as a GRIN (gradient refractive
%   index). This is realized as a set of ellipsoidal iso-indicial surfaces.
%   The particular parameters are taken from Atchison 2006 Vision Research,
%   which were in turn derived from Jones et al 2005 Vision Research and
%   Liou & Brennan 1997 JOSA-A.
%
%   The lens is modeled as being center aligned with the optical axis of
%   the eye, and thus no tilt or shift is needed in the current model. The
%   axial position of the lens center is taken from Atchison 2006.
%
%   The current implementation is a bit of a hodge-podge of different
%   models (Atchison, Navarro, and Jones). It would be nice to implement a
%   more comprehensive model that has an adaptive geometry of the lens
%   interior with accommodative changes, such as:
%
%       Sheil, Conor J., and Alexander V. Goncharov. "Accommodating
%       volume-constant age-dependent optical (AVOCADO) model of the
%       crystalline GRIN lens." Biomedical optics express 7.5 (2016):
%       1985-1999.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   lens                  - Structure.
%


% Because of various imperfections in the model and differences from the
% Navarro paper, it was necessary to "tune" the assigned accomodation
% values so that a requested accomodation state of the emmetropic eye
% results in the expected point of best focus. This is determined in the
% routine 'calcDerivedParams'
accommodationPolyCoef = eye.derivedParams.accommodationPolyCoef;

% The model assumes an 18 year old eye for the Navarro calculations
age = 18;

% Set the D parameter of the model
if ~isempty(eye.meta.navarroD)
    % The D parameter has been hard-coded
    D = eye.meta.navarroD;
else
    if isfield(eye.meta,'accommodationDiopeters')
        accommodationDiopeters = eye.meta.accommodationDiopeters;
    else
        accommodationDiopeters = 0;
    end
    % Convert the requested accommodationDiopeters to the corresponding
    % Navarro D param.
    D = polyval(accommodationPolyCoef,accommodationDiopeters);
end


% Initialize the components of the optical system
lens.S = [];
lens.boundingBox = [];
lens.side = [];
lens.mustIntersect = [];
lens.index = [];
lens.label = {};
lens.plot.color = {};


% Obtain the core and edge refractive indices
nEdge = returnRefractiveIndex( 'lens.edge', eye.meta.spectralDomain, 'age',  age);
nCore = returnRefractiveIndex( 'lens.core', eye.meta.spectralDomain, 'age',  age);

% The position of the point in the lens with the maximal refractive index.
% Positioned so that when the eye is accommodated to focus on a point 200
% mm away (i.e., 5 diopters) the front surface of the lens is at the level
% of the aperture stop of the iris.
lens.center = [-5.6 0 0];

% The thickness of the back and front of the lens, taken from Navarro 2014,
% table 2.
lensThick = 2.93 + 0.0236*age + (0.058 - 0.0005*age)*D;
lensThickBack = 0.6 * lensThick;
lensThickFront = 0.4 * lensThick;
lens.back = lens.center;
lens.back(1) = lens.back(1) - lensThickBack;
lens.front = lens.center;
lens.front(1) = lens.front(1) + lensThickFront;

% The lens has a continuous, gradient refractive index, but is approximated
% by a number of iso-indicial contours, or "shells". The current model is
% only defined for a division of the lens into an odd number of shells.
% Because the lens surfaces are hyperboloids and the gradient shell model
% is ellipsoidal, it can be the case that for some accommodative states,
% the outermost modeled shells do not lie entirely within the lens surface.
% To address this, the current model discards a few of the outermost shells
% so that the gradient set is enclosed entirely within the lens surfaces.
% The number of discarded shells varies by the accommodative state of the
% lens. This has the effect of causing the gradient encountered by the ray
% to be slightly discontinuous close to the outer edges of the lens.
nShells = 21;
endShell = 21;
if D >= 7
    startShell = 3;
end
if D >= 3 && D < 7
    startShell = 4;
end
if D < 4
    startShell = 6;
end


%% Back lens surface
% The back surface of the lens is modeled as a hyperbola. Values taken
% from Navarro 2014. To convert R and Q to radii of a hyperbola:
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
R = 1./( -(1/(5.9 - 0.013*age)) - 0.0043*D  );
Q = -3;
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
radii(1) = abs(b); radii(2:3) = abs(a);

% Build the quadric
S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[lens.center(1)-lensThickBack-radii(1) 0 0]);
boundingBox = [lens.center(1)-lensThickBack lens.center(1) -5 5 -5 5];

% Add to the optical system structure
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; -1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.index = [lens.index; nEdge];
lens.label = [lens.label; {'lens.back'}];
lens.plot.color = [lens.plot.color; {'red'}];

%% Back gradient shells
boundingBox = [lens.center(1)-lensThickBack lens.center(1) -5 5 -5 5];
nLensVals = linspace(nEdge,nCore,nShells*2+1);
for ii = startShell:endShell
    nBackQuadric = [-0.010073731138546 0 0 0.0; 0 -0.002039930555556 0 0; 0 0 -0.002039930555556 0; 0 0 0 1.418000000000000];
    S = nBackQuadric;
    S(end,end)=S(end,end)-nLensVals(ii*2);
    S = quadric.translate(S,[lens.center(1) 0 0]);
    
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
boundingBox = [lens.center(1) lens.center(1)+lensThickFront -5 5 -5 5];
nLensVals = linspace(nEdge,nCore,nShells*2+1);
for ii = endShell:-1:startShell
    nFrontQuadric = [0.022665895061728 0 0 0; 0 0.002039930555556 0 0; 0 0 0.002039930555556 0; 0 0 0 1.371000000000000];
    S = nFrontQuadric;
    S(end,end)=nLensVals(ii*2)-nCore;
    S = quadric.translate(S,[lens.center(1)-0.065277777777778 0 0]);
    
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
% Taken from Navarro 2014.
R = 1/( 1/(12.7-0.058*age) + 0.0077*D  );
Q = -4 - (0.5*D);
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
radii(1) = abs(b);
radii(2:3) = abs(a);

% Build the quadric
S = quadric.scale(quadric.unitTwoSheetHyperboloid, radii);
S = quadric.translate(S,[eye.stop.center(1)+radii(1) 0 0]);
c = quadric.center(S); r = quadric.radii(S);
boundingBox = [lens.center(1) c(1)-r(1) -5 5 -5 5];

% Add to the optical system structure. No refractive index added with this
% surface, as this is the last surface of this set.
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; 1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.label = [lens.label; {'lens.front'}];
lens.plot.color = [lens.plot.color; {'red'}];


end

