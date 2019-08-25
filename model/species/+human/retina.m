function retina = retina(eye)
% Returns the retina sub-field of an eye model structure
%
% Syntax:
%  retina = human.retina( eye )
%
% Description:
%   The retinal surface (i.e., the vitreous chamber / posterior segment) is
%   modeled as an ellipsoidal quadric surface, with adjustments in size for
%   the spherical ametropia of the eye. Atchison 2006 provides the semi-
%   radii of a triaxial ellipsoid model of the retinal surface, with these
%   values varying by spherical ametropia:
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision
%       research 46.14 (2006): 2236-2250.
%
%   The Atchison model is arranged around the visual axis, with the
%   vitreous chamber tilted and shifted. In the current model, the retinal
%   ellipsoid is aligned with and centered on the optical axis.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   retina                - Structure.
%
% Examples:


% Semi-radii of the ellipsoid as a function of ametropia, values taken from
% Atchison 2006.
retinaRadiiEmmetrope = [10.148 11.455 11.365];
retinaRadiiAmetropiaSlope = [-0.163 -0.043 -0.090];
retinaRadii = ...
    retinaRadiiEmmetrope + retinaRadiiAmetropiaSlope.* eye.meta.sphericalAmetropia;

% The current model holds the depth of the anterior chamber constant across
% ametropia, consistent with Atchison.
%
% To position the retina, we need to know the axial length of the eye. If
% this value is not provided, it is derived from the Atchison equation.
if isempty(eye.meta.axialLength)
    axialLength = 23.58 - 0.299*eye.meta.sphericalAmetropia;
else
    axialLength = eye.meta.axialLength;
end
retinaCenter = ...
    [-(axialLength - retinaRadii(1)) 0 0];

% Assemble the components
S = quadric.scale(quadric.unitSphere,retinaRadii);
S = quadric.translate(S,retinaCenter);
retina.S = quadric.matrixToVec(S);
retina.side = -1;
retina.boundingBox = [-30 -7 -15 15 -15 15];
retina.mustIntersect = 1;
retina.label = {'retina'};
retina.plot.color = {[.7,.5,.7]};

end

