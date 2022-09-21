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


% Unless othewise stated, values taken from:
%	Coile, D. C., and L. P. O'Keefe. "Schematic eyes for domestic animals."
%	Ophthalmic and Physiological Optics 8.2 (1988): 215-219.
%
% and
%   Mutti, Donald O., Karla Zadnik, and Christopher J. Murphy. "Naturally
%   occurring vitreous chamber-based myopia in the Labrador retriever."
%   Investigative ophthalmology & visual science 40.7 (1999): 1577-1584.
%
% Values are given for an emmetropic canine eye.


%% Front corneal surface 
% Corneal radii - I don't have a value for asphericity / eccentricity, so
% using the human value
%{
    corneaFrontR = 9.13;
    corneaFrontQ = 0;        
    a = corneaFrontR / ( corneaFrontQ + 1 );
    b = corneaFrontR * sqrt(1/(corneaFrontQ+1)) ;
    radii = [a, b, b]
%}
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


%% Assemble the combined corneal surfaces
cornea.S = [cornea.front.S; cornea.front.S];
cornea.boundingBox = [cornea.front.boundingBox; cornea.front.boundingBox];
cornea.side = [1; 1];
cornea.mustIntersect = [1; 1];
cornea.index = 1.333;
cornea.label = {'cornea.back'; 'cornea.front'};
cornea.plot.color = {[0.5 0.5 0.75]; [0.5 0.5 0.75]};


end