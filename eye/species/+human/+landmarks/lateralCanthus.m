function lateralCanthus = lateralCanthus( eye )
% Returns the lateralCanthus landmark sub-field of an eye model structure
%
% Syntax
%  medialCanthus = human.landmarks.medialCanthus( eye )
%
% Description:
%   Calculates the position of the lateralCanthus.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   lateralCanthus         - Structure with the subfield coords
%
% Examples:
%{
%}



lateralCanthus.coords = [-10 -17 1];

switch eye.meta.eyeLaterality
    case 'Right'
        % No change needed
    case 'Left'
        lateralCanthus.coords(2) = -lateralCanthus.coords(2);
    otherwise
        error('eye laterality not defined')
end


end
