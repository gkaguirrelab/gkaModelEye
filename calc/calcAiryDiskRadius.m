function airyDiskRadius = calcAiryDiskRadius(eye,stopRadius)
% Returns the radius of the Airy Disk for an eye and stop radius
%
% Syntax:
%  airyDiskRadius = calcAiryDiskRadius(eye,stopRadius)
%
% Description
%   The diffraction limit of an optical system can be characterized by the
%   size of the Airy Disk, which is influenced by the wavelength of light
%   and the size of the aperture stop.
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
    diskRadius = nan(size(stopRadius));
    for ss = 1:length(stopRadius)
        diskRadius(ss) = calcAiryDiskRadius(eye,stopRadius(ss));
    end
    plot(stopRadius,diskRadius);
    xlabel('radius aperture stop [mm]');
    ylabel('radius Airy Disk [mm]');
    title('Airy Disk size by stop radius');
%}


arguments
    eye (1,1) {isstruct}
    stopRadius (1,1) {mustBeNumeric}
end

% The wavelength for the calculation, in mm
lambdaMm = eye.index.wavelength / 1e6;

% The angular width of the first minimum of the Airy Disk on the retina,
% as seen from the aperture, in degrees
angularDiam = rad2deg(1.22 * lambdaMm / (stopRadius*2));

% We obtain the distance of the aperture stop from the vertex of the retina
d = eye.stop.center(1) - eye.landmarks.vertex.coords(1);

% Convert the angular diameter to radius, and then mm
airyDiskRadius = tand(angularDiam/2)*d;

end
