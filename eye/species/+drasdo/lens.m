function lens = lens( eye )
% Returns the lens sub-field of an eye model structure
%
% Syntax:
%  lens = human.lens( eye )
%
% Description:

% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   lens                  - Structure.
%

%% Lens front surface
radii = [10 10 10];

% Create the quadric
Sfront = quadric.scale(quadric.unitSphere,radii);

% Shift the quadric so that the apex of the front surface is at -3.6 mm
anteriorLensPosition = -3.6;
Sfront = quadric.translate(Sfront,[-radii(1)+anteriorLensPosition 0 0]);

% Set the bounding box
bbfront = [-5 anteriorLensPosition -5 5 -5 5];


%% Lens back surface
radii = [6 6 6];

% Create the quadric
Sback = quadric.scale(quadric.unitSphere,radii);

% Shift the quadric so that the posterior apex of the back surface is at
% -7.375 mm. This location was not provided in the Drasdo & Fowler 1974
% paper, but was obtained by measurement of the figure by Gionvanni
% Montesano, and reported here:
%
%   Montesano, Giovanni, et al. "Revisiting the drasdo model: implications
%   for structure-function analysis of the macular region." Translational
%   Vision Science & Technology 9.10 (2020): 15-15.
% 
posteriorLensPosition = -7.375;
Sback = quadric.translate(Sback,[radii(1)+posteriorLensPosition 0 0]);

% Set the bounding box
bbback = [posteriorLensPosition -5 -5 5 -5 5];


%% Add some locations on the front, middle, and back of the lense
% These are used as starting points in ray tracing routines
lens.back = [posteriorLensPosition 0 0];
lens.center = [-5 0 0];
lens.front = [anteriorLensPosition 0 0];


%% Assemble the system
lens.S = [quadric.matrixToVec(Sback); quadric.matrixToVec(Sfront)];
lens.boundingBox = [bbback; bbfront];
lens.side = [-1; 1];
lens.mustIntersect = [1; 1];
lens.index = [1.430];
lens.label = [{'lens.back'};{'lens.front'}];
lens.plot.color = [{[0.75 0.75 0.75 0.75]};{[0.75 0.75 0.75 0.75]}];
lens.meta.navarroD = nan;

end

