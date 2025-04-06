function stop = stop( eye )
% Returns the stop sub-field of an eye model structure
%
% Syntax:
%  stop = human.stop( eye )
%
% Description:
%   The stop is an aperture in the iris, centered on the longitudinal axis.
%   The stop is modeled as an ellipse, with the eccentricity and theta
%   varying with dilation.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   stop                  - Structure.
%

% Center the stop on the longitudinal axis within the iris
stop.center = [eye.iris.center(1) 0 0];

% The aperture stop of the eye is elliptical. The parameters of this 
% are determined in the routine 'calcDerivedParams', and passed here in the
% eye.derivedParams field.
stop.eccenParams = eye.derivedParams.stopEccenParams;

% Specify the params and equation that defines the stop ellipse.
% This can be invoked as a function using str2func.
stop.eccenFcnString = sprintf('@(x) (tanh((x+%f).*%f)+%f)*%f',stop.eccenParams(1),stop.eccenParams(2),stop.eccenParams(3),stop.eccenParams(4));

% The theta values of the stop ellipse for eccentricities less
% than, and greater than, zero.
switch eye.meta.eyeLaterality
    case 'Right'
        stop.thetas = [0  3/7*pi];
    case 'Left'
        stop.thetas = [0  4/7*pi];
end

end

