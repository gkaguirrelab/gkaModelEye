function stop = stop( eye )
% Returns the stop sub-field of an eye model structure
%
% Syntax:
%  stop = human.stop( eye )
%
% Description:
%   The stop is an aperture in the iris, centered on the optical axis. The
%   stop is modeled as an ellipse, with the eccentricity and theta varying
%   with dilation.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   stop                  - Structure.
%

% Center the stop on the optical axis within the iris
stop.center = [-3.6 0 0];

% Make the stop circular.
stop.eccenParams = [0 0 0 0];

% Specify the params and equation that defines the stop ellipse.
% This can be invoked as a function using str2func.
stop.eccenFcnString = sprintf('@(x) 0');

% The theta values of the stop ellipse for eccentricities less
% than, and greater than, zero.
switch eye.meta.eyeLaterality
    case 'Right'
        stop.thetas = [0  0];
    case 'Left'
        stop.thetas = [0  0];
end

end

