function mmPerDeg = calcMmRetinaPerDeg(eye,fieldLocation,deltaDegEuclidean,cameraMedium)
% Returns mm/degree of visual angle at the given visual field location
%
% Syntax:
%  mmPerDeg = calcMmRetinaPerDeg(eye,fieldLocation,deltaDegEuclidean,cameraMedium)
%
% Description
%   Given an eye structure and the location of a point in visual space
%   (w.r.t. the optical axis of the eye), returns the number of mm of
%   retina per degree of visual angle at the corresponding retinal location
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   degField              - 2x1 vector that specifies the visual field
%                           position of a point in degree angles along the
%                           horizontal and vertical meridian with respect
%                           to the optical axis of the eye.
%   deltaDegEuclidean     - Scalar. The measurement is made for a small
%                           displacement of visual angle specified by this
%                           variable.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   mmPerDeg              - Scalar. The mm of retina per degree of visual
%                           angle at the specified visual field location.
%
% Examples:
%{
    % mm of retina / deg at the fovea
    eye = modelEyeParameters('calcLandmarkFovea',true);
    fieldLocation = eye.landmarks.fovea.degField(1:2);
    mmPerDeg = calcMmRetinaPerDeg(eye,fieldLocation);
    fprintf('%2.3f retinal mm per deg visual field at the fovea in the emmetropic eye.\n',mmPerDeg);
%}

% Handle incomplete inputs
if nargin==0
    error('Need to specify an eye structure');
end

if nargin==1
    fieldLocation = [0 0];
    deltaDegEuclidean = 1e-2;
    cameraMedium = 'air';
end

if nargin==2
    deltaDegEuclidean = 1e-2;
    cameraMedium = 'air';
end

if nargin==3
    cameraMedium = 'air';
end

deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];

[~,X0] = calcRetinaFieldPoint( eye, fieldLocation-deltaAngles./2, cameraMedium);
[~,X1] = calcRetinaFieldPoint( eye, fieldLocation+deltaAngles./2, cameraMedium);

mmPerDeg = norm(X0-X1) / norm(deltaAngles);


end