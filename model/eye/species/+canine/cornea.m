function cornea = cornea( eye )
% Returns the cornea sub-field of an eye model structure
%
% Syntax:
%  cornea = human.cornea( eye )
%
% Description:
% 
% Unless othewise stated, values taken from:
%   Mutti, Donald O., Karla Zadnik, and Christopher J. Murphy. "Naturally
%   occurring vitreous chamber-based myopia in the Labrador retriever."
%   Investigative ophthalmology & visual science 40.7 (1999): 1577-1584.
%
% Values are given for an emmetropic canine eye.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   cornea                - Structure.
%


%% Front corneal surface 
radii = [7.1 7.1 7.1];

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
cornea.front.boundingBox=[-4.84 X(1) -9 9 -9 9];
cornea.front.center=[-radii(1) 0 0];

%% Back corneal surface 
% Model this as the same as the front. Given that the index of refraction
% of the cornea and the aqueous is the same in this reduced eye, this has
% no effect upon the optical properties
radii = [7.1 7.1 7.1];

% Create the quadric
S = quadric.scale(quadric.unitSphere,radii);

corneaThickness = 0.01;

% We set the center of the cornea front surface ellipsoid so that the axial
% apex (prior to rotation) is at position [0, 0, 0]
S = quadric.translate(S,[-radii(1)-corneaThickness 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store these values
cornea.front.S = quadric.matrixToVec(S);
cornea.front.side = 1;
cornea.front.boundingBox=[-4.84 X(1) -9 9 -9 9];
cornea.front.center=[-radii(1) 0 0];


%% Assemble the combined corneal surfaces
cornea.S = [cornea.front.S; cornea.front.S];
cornea.boundingBox = [cornea.front.boundingBox; cornea.front.boundingBox];
cornea.side = [1; 1];
cornea.mustIntersect = [1; 1];
cornea.index = 1.333;
cornea.label = {'cornea.back'; 'cornea.front'};
cornea.plot.color = {[0.5 0.5 0.75]; [0.5 0.5 0.75]};


end