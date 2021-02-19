function mmPerDeg = calcMmRetinaPerDeg(eye,fieldOrigin,deltaDegEuclidean,rayOriginDistance,cameraMedium)
% Returns mm/degree of visual angle at the given visual field location
%
% Syntax:
%  mmPerDeg = calcMmRetinaPerDeg(eye,fieldOrigin,deltaDegEuclidean,cameraMedium)
%
% Description
%   Given an eye structure and the location of a point in visual space
%   (w.r.t. the optical axis of the eye), returns the number of mm of
%   retina per degree of visual angle at the corresponding retinal location
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   fieldOrigin           - 1x2 or 2x1 vector that provides the coordinates
%                           in degrees of visual angle of the origin of the
%                           nodal ray.
%   deltaDegEuclidean     - Scalar. The measurement is made for a small
%                           displacement of visual angle specified by this
%                           variable.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex. Assumed to be
%                           500 mm if not defined.
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
    eye = modelEyeParameters();
    fieldOrigin = eye.landmarks.fovea.degField(1:2);
    mmPerDeg = calcMmRetinaPerDeg(eye,fieldOrigin);
    fprintf('%2.3f retinal mm per deg visual field at the fovea in the emmetropic eye.\n',mmPerDeg);
%}

% Handle incomplete inputs
if nargin==0
    error('Need to specify an eye structure');
end

if nargin==1
    fieldOrigin = [0 0];
    deltaDegEuclidean = 1e-2;
    rayOriginDistance = 1000;
    cameraMedium = 'air';
end

if nargin==2
    deltaDegEuclidean = 1e-2;
    rayOriginDistance = 1000;
    cameraMedium = 'air';
end

if nargin==3
    rayOriginDistance = 1000;
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
end

% Set up the jitter in angles around the specified field location
deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];

rayPath0 = calcNodalRayFromField(eye,fieldOrigin-deltaAngles./2,rayOriginDistance,cameraMedium);
rayPath1 = calcNodalRayFromField(eye,fieldOrigin+deltaAngles./2,rayOriginDistance,cameraMedium);

mmPerDeg = norm(rayPath0(:,end)-rayPath1(:,end)) / norm(deltaAngles);


end