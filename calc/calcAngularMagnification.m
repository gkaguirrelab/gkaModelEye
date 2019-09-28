function M = calcAngularMagnification(eye, varargin)
% Calcuates the percept angular magnification produced by artificial lenses 
%
% Syntax:
%  magnification = calcAngularMagnification(eye)
%
% Description
%   Calculates the angular magnification as perceived by the eye produced
%   the addition of artificial lenses (contacts or spectacles). Values
%   greater than 1 indicate magnification, values less than one reflect
%   minification.
%
%   
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%
% Optional key/value pairs:
%   All varargin are passed on to assembleOpticalSystem. For this routine,
%   typical options are:
%
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, or 1x3 vector, with values for the
%                           lens refraction in diopters, (optionally) the
%                           index of refraction of the lens material, and
%                           (optinally) the vertex distance in mm. If left
%                           empty, no spectacle is added to the model.
%
% Outputs:
%   M                     - Scalar. The magnification produced by the
%                           system.
%
% Examples:
%{
    eye = modelEyeParameters;
    calcAngularMagnification(eye,'spectacleLens',-4)
%}

% A point on the retina, just offset from the vertex
X = eye.landmarks.vertex.coords';
X(2) = 0.1;

% The optical center of the eye
P = calcOpticalCenter(eye);

% A ray departing the retina, aimed at the optical center
R = quadric.normalizeRay([X P-X]);

% Optical systems with and without the artificial lens
opticalSystemNoLens = assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[]);
opticalSystemWithLens = assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[],varargin{:});

% Ray trace 
outputRayNoLens = rayTraceQuadrics(R,opticalSystemNoLens);
outputRayWithLens = rayTraceQuadrics(R,opticalSystemWithLens);

% Find the ratio of the angle w.r.t. the optical axis for the without and
% with lens case. This is the magnification.
M = quadric.rayToAngles(outputRayNoLens) / quadric.rayToAngles(outputRayWithLens);

end

