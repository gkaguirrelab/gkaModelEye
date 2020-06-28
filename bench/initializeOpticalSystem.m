function opticalSystem = initializeOpticalSystem(refractiveIndex)
% Returns the first row of an optical system
%
% Syntax:
%  opticalSystem = initializeOpticalSystem(refractiveIndex)
%
% Description:
%   The optical system matrix begins with a first row that is all nans,
%   with the exception of the last value which is the refractive index of
%   the medium in which the ray arises. This simple routine returns this
%   first row, populated with a specified refractive index.
%
% Inptuts:
%   refractiveIndex       - Scalar. The refractive index of the medium
%
% Outputs:
%   opticalSystem         - 1x19 matrix. This is the first row of a valid
%                           optical system.
%

% If no refractive index value was passed, assumed the a value of 1 for
% air.
if nargin == 0
    refractiveIndex = 1;
end

% The opticalSystem starts with a row of nans...
opticalSystem = nan(1,19);

% ...except for the refractive index value in the last position.
opticalSystem(1,19) = refractiveIndex;

end