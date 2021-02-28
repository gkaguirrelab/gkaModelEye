function [incidentNode,emergentNode] = nodes( eye )
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
[iN,eN] = calcNodes(eye);

% Create the structure to return
incidentNode.coords = iN';
emergentNode.coords = eN';

end
