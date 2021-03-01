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

% By convention, the position of the nodes are calculated and stored for
% the eye when the navarroD parameter of the lens is set to zero.
eye.meta.accommodation = [];
eye.meta.navarroD = 0;
eye.lens = human.lens(eye);

% Obtain the cartesian coordinates of the nodes
[iN,eN] = calcNodes(eye);

% Create the structure to return
incidentNode.coords = iN';
incidentNode.meta = 'Calculated for lens with navarroD = 0';
emergentNode.coords = eN';
emergentNode.meta = 'Calculated for lens with navarroD = 0';

end
