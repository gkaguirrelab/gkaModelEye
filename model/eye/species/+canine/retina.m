function retina = retina(eye)
% Returns the retina sub-field of an eye model structure
%
% Syntax:
%  retina = human.retina( eye )
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
%   retina                - Structure.
%
% Examples:



% An essentially spherical retina. I add a tiny bit of asphericity to allow
% the ellipsoidal coordinate system to work with this model.
radii = [ 9.7 9.701 9.702];

% Create the quadric
S = quadric.scale(quadric.unitSphere,radii);

% Shift so that the posterior apex of the sphere is 21.5 mm posterior to
% the corneal apex
S = quadric.translate(S,[radii(1)-21.5 0 0]);

% Assemble the system
retina.S = quadric.matrixToVec(S);
retina.side = -1;
retina.boundingBox = [-21.5 -4.84 -12 12 -12 12];
retina.mustIntersect = 1;
retina.label = {'retina'};
retina.plot.color = {[0.75, 0.75, 0.75]};

end

