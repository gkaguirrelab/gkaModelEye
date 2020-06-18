function medialCanthus = medialCanthus( eye )
% Returns the medial canthus landmark sub-field of an eye model structure
%
% Syntax
%  medialCanthus = human.landmarks.medialCanthus( eye )
%
% Description:
%   Calculates the position of the medial canthus.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   medialCanthus         - Structure with the subfield coords
%
% Examples:
%{
%}



medialCanthus.coords = [-3 12.5 -1];

switch eye.meta.eyeLaterality
    case 'Right'
        % No change needed
    case 'Left'
        medialCanthus.coords(2) = -medialCanthus.coords(2);
    otherwise
        error('eye laterality not defined')
end


end
