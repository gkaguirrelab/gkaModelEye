function nodes = nodes( eye )
% Returns the nodal points landmark sub-field of an eye model structure
%
% Syntax
%  nodalPoints = human.landmarks.nodalPoints( eye )
%
% Description:
%   Calculates the approximate incident and emergent nodal point of the
%   eye.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   nodalPoints           - Structure with the subfield coords
%
% Examples:
%{
%}

% Obtain the cartesian coordinates of the nodes
[incidentNode,emergentNode] = calcNodes(eye);

% Create the structure to return
nodes.incident.coords = incidentNode';
nodes.emergent.coords = emergentNode';


end
