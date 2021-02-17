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
Sfront = quadric.translate(Sfront,[-radii(1)-3.6 0 0]);

% Set the bounding box
bbfront = [-5 -3.6 -5 5 -5 5];


%% Lens back surface
radii = [6 6 6];

% Create the quadric
Sback = quadric.scale(quadric.unitSphere,radii);

% Shift the quadric so that the posterior apex of the back surface is at
% -7.375 mm
Sback = quadric.translate(Sback,[radii(1)-7.375 0 0]);

% Set the bounding box
bbback = [-7.4 -5 -5 5 -5 5];


%% Add some locations on the front, middle, and back of the lense
% These are used as starting points in ray tracing routines
lens.back = [-7.4 0 0];
lens.center = [-5 0 0];
lens.front = [-3.6 0 0];


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

