function n = returnRefractiveIndex( material, lambda, varargin )
% Refractice index for a specified material at a specified wavelength
%
% Syntax:
%  n = returnRefractiveIndex( material, wavelength )
%
% Description:
%   The ray tracing model requires the refractive index of several
%   biological and optical materials. The index of refraction varies by
%   wavelength of light. This routine returns the index for a specified
%   material for a specified imaging domain (visual or near infrared).
%
%   The refractive index for a material as a function of wavelength is
%   given by the Cauchy equation, with the values for biological materials
%   taken from:
%
%       Navarro, Rafael. "Adaptive model of the aging emmetropic eye and
%       its changes with accommodation." Journal of vision 14.13 (2014):
%       21-21.
%
%
% Inputs:
%   material              - Char vector.
%   lambda                - Scalar or char vector. Identifies the 
%                           wavelength (in nm) for which the refractive
%                           index should be calculated. If char, valid
%                           values are {'VIS','NIR'}.
%
% Optional key/value pairs:
%  'age'                  - Scalar, age in years. Defaults to 18. Used in
%                           the calculation of the parameters of the lens
%                           core.
%
% Outputs:
%   n                     - The requested index of refraction.
%
% Examples:
%{
    n = returnRefractiveIndex( 'cornea', 'visible' )
%}

%% input parser
p = inputParser;

% Required
p.addRequired('material',@ischar);
p.addRequired('lambda',@(x)(ischar(x) || isscalar(x)));

% Optional
p.addParameter('age',18,@isscalar);

% parse
p.parse(material,lambda, varargin{:});



% Check or assign lambda
if ischar(p.Results.lambda)
    switch p.Results.lambda
        case {'VIS','vis','Vis','visible'}
            % Set lambda to 555 nm, the peak of the photopic luminosity
            % function
            lambda = 555;
        case {'NIR','nir','Nir','near infrared','near infra-red'}
            % Set lambda to 950 nm, which is the center frequency of light
            % emitted by the LED of many active IR cameras
            lambda = 950;
        otherwise
            error(['I do not have values for the spectral domain of ' spectralDomain]);
    end
end

% Assign the coefficients of the Cauchy equation
switch material
    case 'vacuum'
        c = [1 0 0 0];
    case 'air'
        c = [1 0 0 0];
    case 'water'
        % Bashkatov, Alexey N., and Elina A. Genina. "Water refractive
        % index in dependence on temperature and wavelength: a simple
        % approximation." Saratov Fall Meeting 2002: Optical Technologies
        % in Biophysics and Medicine IV. Vol. 5068. International Society
        % for Optics and Photonics, 2003.
        c = [1.3176, 5.51547658e3, -2.5756e8 9.47474];
    case 'tears'
        % Patel, Sudi, Karen E. Boyd, and Janet Burns. "Age, stability of
        % the precorneal tear film and the refractive index of tears."
        % Contact Lens and Anterior Eye 23.2 (2000): 44-47.
        c = [1.323757, 5.560240e3, -5.817391e8, 5.036810e13];
    case 'vitreous'
        c = [1.323757, 5.560240e3, -5.817391e8, 5.036810e13];
    case 'lens.core'
        A = 1.40965 + (3.55e-4 * p.Results.age) - (7.5e-6 * p.Results.age);
        c = [A, 6.521218e3, -6.11066e8, 5.908191e13];      
    case 'lens.edge'
        c = [1.356086, 6.428455e3, -6.023738e8, 5.824149e13];
    case 'aqueous'
        c = [1.323031, 6.070796e3, -7.062305e8, 6.147861e13];
    case 'cornea'
        c = [1.362994, 6.009687e3, -6.760760e8, 5.908450e13];
    case 'hydrogel'
        % Childs, Andre, et al. "Fabricating customized hydrogel contact
        % lens." Scientific reports 6 (2016): 34905.
        ns = [1.430 1.420];
    case 'optorez'
        % Nikolov, Ivan D., and Christo D. Ivanov. "Optical plastic
        % refractive measurements in the visible and the near-infrared
        % regions." Applied Optics 39.13 (2000): 2067-2070.
        % Measured at 594 and 890 nm
        ns = [1.5089 1.5004];
    case 'cr-39'
        % Traynor, Nathan BJ, et al. "CR-39 (PADC) Reflection and
        % Transmission of Light in the Ultraviolet-Near-Infrared (UV-NIR)
        % Range." Applied spectroscopy (2017): 0003702817745071.
        % NIR value estimated from Figure 5b.
        ns = [1.51 1.50];
    case 'polycarbonate'
        % Nikolov, Ivan D., and Christo D. Ivanov. "Optical plastic
        % refractive index measurements for NIR region." 18th Congress of
        % the International Commission for Optics. Vol. 3749. International
        % Society for Optics and Photonics, 1999.
        ns = [1.5852 1.5614];
    otherwise
        error(['I do not know the index of refraction for ' material]);
end

% Cauchy?s (1836) equation (cited in Atchison & Smith 2005)
n = c(1) + c(2)/lambda^2 + c(3)/lambda^4 + c(4)/lambda^6;



end

