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

% By convention, the position of the nodes are calculated and stored with
% lens of the eye set to a state that would provide accommodation at
% infinity for the default, emmetropic eye.
eye.meta.accommodation = [];
eye.meta.navarroD = eye.derivedParams.navarroDAtInfinity;
eye.lens = human.lens(eye);

% Obtain the cartesian coordinates of the nodes
[iN,eN] = calcNodes(eye);

% Create the structure to return
incidentNode.coords = iN';
incidentNode.meta = sprintf('Calculated for lens with navarroD = %2.2f',eye.derivedParams.navarroDAtInfinity);
emergentNode.coords = eN';
emergentNode.meta = sprintf('Calculated for lens with navarroD = %2.2f',eye.derivedParams.navarroDAtInfinity);

end
