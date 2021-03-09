function airyDiskRadius = calcAiryDiskRadius(eye,stopRadius)
% Returns the radius of the Airy Disk for an eye and stop radius
%
% Syntax:
%  airyDiskRadius = calcAiryDiskRadius(eye,stopRadius)
%
% Description
%   Foo
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   stopRadius            - Scalar. The size in mm of the aperture stop.
%
% Outputs:
%   airyDiskRadius        - Scalar. Radius (in mm) of the Airy Disk.
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters('spectralDomain',550,'accommodation',1);
    % Create a plot of Airy Disk diameters vs. stop size
    stopRadius = 0.5:0.5:3;
    diskDiam = [];
    for ss = stopRadius
        diskDiam(end+1) = calcAiryDiskRadius(eye,ss);
    end
    plot(stopRadius,diskDiam);
%}


arguments
    eye (1,1) {isstruct}
    stopRadius (1,1) {mustBeNumeric} = 1
end

% Obtain the optical power for this eye
%% CALCULATE FOR FOCAL LENGTH IN AIR OR IN VITREOUS?
[diopters,focalPoint] = calcOpticalPower(eye);

% Obtain the radius of the entrance pupil for this eye
pupilRadius = stopRadius;

% The f-number of the eye
f = norm(1000/diopters)/(pupilRadius*2);

% The wavelength for the calculation, in mm
lambdaMm = eye.index.wavelength / 1e6;

% The airy disk diameter
airyDiskRadius = 2.44 * lambdaMm * f;


end
