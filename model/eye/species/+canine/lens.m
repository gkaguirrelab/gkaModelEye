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
radii = [4.5 4.5 4.5];

% Create the quadric
Sfront = quadric.scale(quadric.unitSphere, radii);

anteriorLensPosition = -4.27;
Sfront = quadric.translate(Sfront,[-radii(1)+anteriorLensPosition 0 0]);

% Set the bounding box
bbfront = [-9 anteriorLensPosition -5 5 -5 5];


%% Lens back surface
%{
    R = -6.520;
    Q = -2;
    a = R * sqrt(abs( 1 / (Q - 1 ) )) * sign(Q);
    b = R / (Q - 1 );
    radii = [a b b]
%}
radii = [5.2 5.2 5.2];

% Create the quadric
Sback = quadric.scale(quadric.unitSphere, radii);

posteriorLensPosition = -11.54;
Sback = quadric.translate(Sback,[radii(1)+posteriorLensPosition 0 0]);

% Set the bounding box
bbback = [posteriorLensPosition -9 -5 5 -5 5];


%% Assemble the system
lens.S = [quadric.matrixToVec(Sback); quadric.matrixToVec(Sfront)];
lens.boundingBox = [bbback; bbfront];
lens.side = [-1; 1];
lens.mustIntersect = [1; 1];
lens.index = [1.5361];
lens.label = [{'lens.back'};{'lens.front'}];
lens.plot.color = [{[0.75 0.75 0.75 0.75]};{[0.75 0.75 0.75 0.75]}];
lens.meta.navarroD = nan;

end

