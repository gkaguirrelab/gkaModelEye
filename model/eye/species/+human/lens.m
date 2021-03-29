function lens = lens( eye )
% Returns the lens sub-field of an eye model structure
%
% Syntax:
%  lens = human.lens( eye )
%
% Description:
%   A model of the crystalline lens is generated, expressed as a set of
%   quadric surfaces. The anterior and posterior surfaces of the lens are
%   modeled as one-half of a two-sheeted hyperboloid.
%
%   The lens model begins with published values, which are subsequently
%   adjusted to allow for relatively equal accommodation across the retinal
%   surface in the emmetropic eye.
%
%   The starting point for the radii of the surfaces, and their dependence
%   upon age and the accommodative state of the eye, are taken from:
%
%       Navarro, Rafael. "Adaptive model of the aging emmetropic eye and
%       its changes with accommodation." Journal of vision 14.13 (2014):
%       21-21.
%
%   The total lens thickness, and thickness of the anterior and posterior
%   portions of the lens, are also taken from Navarro 2014.
%
%   The interior of the lens is modeled as a GRIN (gradient refractive
%   index). This is realized as a set of iso-indicial surfaces. The
%   particular parameters are taken from Atchison 2006 Vision Research,
%   which were in turn derived from Jones et al 2005 Vision Research and
%   Liou & Brennan 1997 JOSA-A.
%
%   The lens is modeled as being center aligned with the optical axis of
%   the eye, and thus no tilt or shift is needed in the current model.
%
%   The axial position of the lens center was selected so that the anterior
%   surface of the lens is at the leve of the aperture stop when an
%   emmetropic eye is accommodated to 5 D.
%
%   The values for the radius of curvature and asphericity of the posterior
%   surface of the lens were modified so that the image surface (±40°)
%   produced by an emmetropic eye accommodated at infinity closely
%   approximates the retinal surface. SEE: d004_imageSurface.mlx
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   lens                  - Structure.
%


%% Set accommodation
% The accommodative state of the eye can be specified in two different
% ways. An explicit value for the "navarroD" parameter may be set, or the
% desired accommodative state of the eye may be specified. If neither is
% set, then navarroD is set to the value in the derived parameters that
% provides for accommodation at infinity for the default (emmetropic) eye.
% If the navarroD is set, then this value is used. If the accommodation is
% set, then a search is performed in an attempt to find the navarroD
% parameter that provides the desired accommodation value for this eye
% along the longitudinal axis. If both values are set, then the navarroD
% parameter is used as the initial guess for the search to produce the
% desired accommodation.

% Set the navarroD to empty, or the passed value
D = [];
if isfield(eye.meta,'navarroD')
    if ~isempty(eye.meta.navarroD)
        D = eye.meta.navarroD;
    end
end

% Set the desired accommodation to empty, or the passed value
A = [];
if isfield(eye.meta,'accommodation')
    if ~isempty(eye.meta.accommodation)
        A = eye.meta.accommodation;
    end
end

% If neither A nor D was provided, set D to the derived parameter value
if isempty(D) && isempty(A)
    D = eye.derivedParams.navarroDAtInfinity;
end

% If an A value is specified, search to set a D value that provides the
% desired A accommodation.
if ~isempty(A)
    
    % We need to set an initial value for the lens to start the search. To
    % do so, we create a "bootstrap" version of the eye that is initialized
    % with a value for navarroD that provides for resting accommodation. In
    % the bootstrap eye, we also need to empty the "accommodation" field to
    % avoid infinite recursion.
    bootstrapEye = eye;
    bootstrapEye.meta.accommodation = [];
    if isempty(D)
        bootstrapEye.meta.navarroD = eye.derivedParams.navarroDAtInfinity;
    else
        bootstrapEye.meta.navarroD = D;
    end
    bootstrapEye.lens = human.lens(bootstrapEye);
    
    % Now perform the search to find the navarroD parameter that provides
    % the desired accommodation along the longitudinal axis.
    D = calcNavarroD(bootstrapEye,A);
    
    % Done with the bootstrap
    clear bootstrapEye
end


%% Prepare lens parameters

% Obtain the age of the modeled subject
age = eye.meta.ageYears;

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

% The refractive index varies as a function of relative position between
% the edge (1) and core (0). This anonymous function follows Atchison 2006,
% and specifies the relative position of a gradient shell as a function of
% the refractive index of that shell.
nPosition = @(n) (nCore - n).^(1/2)./(nCore - nEdge).^(1/2);

% The position of the point in the lens with the maximal refractive index.
% Positioned so that when the eye is accommodated to focus on a point 200
% mm away (i.e., 5 diopters) the front surface of the lens is at the level
% of the aperture stop of the iris.
lens.center = [-5.2 0 0];

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
numShells = 21;


%% Back lens surface
% The back surface of the lens is modeled as a hyperbola. Values taken
% from Navarro 2014. To convert R and Q to radii of a hyperbola:
%   R = b^2/a
%	Q = (a^2 / b^2) + 1
% Therefore, given R and Q, we can obtain a and b, which correspond
% to the radii of the hyperbola model, with a corresponding to the
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
% Rconst is 5.9 if derived from the Navarro 2014 parameters, but is
% modified here to a value of 4.0 to provide appropriate peripheral retinal
% accommodation in the emmetropic eye.
Rconst = 4.0;
R = 1./( -(1/(Rconst - 0.013*age)) - 0.0043*D  );

% The asphericity (Q) value is given as -3 in Navarro 2014, but is adjusted
% here to be -2.5 to provide appropriate peripheral retinal accommodation
% in the emmetropic eye.
Q = -2.5;
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
lensBackRadii(1) = abs(b); lensBackRadii(2:3) = abs(a);

% Build the quadric. We save the back lens surface primitive for use below
% in the gradient refractive index
backQuadricPrimitive = quadric.scale(quadric.unitTwoSheetHyperboloid, lensBackRadii./lensBackRadii(1));
S = quadric.scale(backQuadricPrimitive,lensBackRadii(1));
S = quadric.translate(S,[lens.back(1)-lensBackRadii(1) 0 0]);

% A bounding box for the back lens surface
boundingBox = [lens.back(1) lens.center(1) -5 5 -5 5];

% Add to the optical system structure
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; -1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.index = [lens.index; nEdge];
lens.label = [lens.label; {'lens.back'}];
lens.plot.color = [lens.plot.color; {[0.75 0.75 0.75 0.75]}];


%% Back gradient shells
% The shells are scaled and then shifted to be equally spaced across
% the thickness of the lens.

% Create a set of equally spaced refractive indices
nLensVals = linspace(nEdge,nCore,numShells*2+1);

for ii = 1:numShells
    
    % Find the position of this shell within the back half of the lens,
    % where position 1 is the edge and 0 is the core
    proportion = nPosition(nLensVals(ii));
    
    % The radius of the quadric shell is equal to proportionally scaled
    % form of the lens back radius
    radius = lensBackRadii(1) * proportion;
    
    % Scale the quadric surface to this radius size    
    S = quadric.scale(backQuadricPrimitive,radius);

    % Position the surface 
    S = quadric.translate(S,[lens.center(1)-(proportion*lensThickBack)-radius 0 0]);
    
    % Create a bounding box
    boundingBox = [lens.center(1)-(proportion*lensThickBack) lens.center(1) -5 5 -5 5];
    
    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; -1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii*2+1)];
    lens.label = [lens.label; {sprintf('lens.back.shell_n=%0.3f',nLensVals(ii*2))}];
    lens.plot.color = [lens.plot.color; [0.75 0.75+0.25*(ii/numShells) 0.75 0.1]];
end
% Force the index of the center shell to be equal to the lens core
% refractive index. The value can be something other than this when
% endShell shell is not equal to the nShells.
lens.index(end) = nCore;


%% Front lens surface
% This is calculated now so that the form is available for the
% construction of the front gradient shells, but not stored in the optical
% system until after the lens gradient is finished

% Taken from Navarro 2014.
R = 1/( 1/(12.7-0.058*age) + 0.0077*D  );
Q = -4 - (0.5*D);
a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
b = R / (Q - 1 );
lensFrontRadii(1) = abs(b);
lensFrontRadii(2:3) = abs(a);

% Build the quadric. We save the back lens surface primitive for use below
% in the gradient refractive index
frontQuadricPrimitive = quadric.scale(quadric.unitTwoSheetHyperboloid, lensFrontRadii./lensFrontRadii(1));
S = quadric.scale(frontQuadricPrimitive,lensFrontRadii(1));
S = quadric.translate(S,[lens.front(1)-lensFrontRadii(1) 0 0]);


%% Front gradient shells
% The shells are scaled and then shifted to be equally spaced across
% the thickness of the lens.

% Create a set of equally spaced refractive indices
nLensVals = linspace(nEdge,nCore,numShells*2+1);

for ii = numShells:-1:1
    
    % Find the position of this shell within the front half of the lens,
    % where position 1 is the edge and 0 is the core
    proportion = nPosition(nLensVals(ii));
    
    % The radius of the quadric shell is equal to proportionally scaled
    % form of the lens front radius
    radius = lensFrontRadii(1) * proportion;
    
    % Scale the quadric surface to this radius size    
    S = quadric.scale(frontQuadricPrimitive,radius);

    % Position the surface 
    S = quadric.translate(S,[lens.center(1)+(proportion*lensThickFront)+radius 0 0]);

    % Create a bounding box
    boundingBox = [lens.center(1) lens.center(1)+(proportion*lensThickFront) -5 5 -5 5];

    % Add this shell to the optical system structure
    lens.S = [lens.S; quadric.matrixToVec(S)];
    lens.boundingBox = [lens.boundingBox; boundingBox];
    lens.side = [lens.side; 1];
    lens.mustIntersect = [lens.mustIntersect; 0];
    lens.index = [lens.index; nLensVals(ii*2-1)];
    lens.label = [lens.label; {sprintf('lens.front.shell_n=%0.3f',nLensVals(ii*2))}];
    lens.plot.color = [lens.plot.color; [0.75 0.75+0.25*(ii/numShells) 0.75 0.1]];
end

% Force the index of the space between the last front shell and the front
% lens surface to be equal to the lens edge refractive index. The value can
% be something other than this when the start shell is something other than
% one.
lens.index(end) = nEdge;

% Add the front surface to the optical system structure. No refractive index added with this
% surface, as this is the last surface of this set.
boundingBox = [lens.center(1) lens.front(1) -5 5 -5 5];
lens.S = [lens.S; quadric.matrixToVec(S)];
lens.boundingBox = [lens.boundingBox; boundingBox];
lens.side = [lens.side; 1];
lens.mustIntersect = [lens.mustIntersect; 1];
lens.label = [lens.label; {'lens.front'}];
lens.plot.color = [lens.plot.color; {[0.75 0.75 0.75 0.75]}];
lens.meta.navarroD = D;

end

