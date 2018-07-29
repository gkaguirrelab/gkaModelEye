function retina = retina(eye)
% Returns the retina sub-field structure of an eye model structure
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
% The Atchison model is arranged around the visual axis, with the vitreous
% chamber tilted and shifted. In the current model, the retinal ellipsoid
% is aligned with and centered on the optical axis.
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
retinaRadiiEmetrope = [10.148 11.455 11.365];
retinaRadiiAmetropiaSlope = [-0.163 -0.043 -0.090];
retinaRadii = ...
    retinaRadiiEmetrope + retinaRadiiAmetropiaSlope.* eye.meta.sphericalAmetropia;

% The current model holds the depth of the anterior chamber constant across
% ametropia, consistent with Atchison, although see:
%
%   Hosny, Mohamed, et al. "Relationship between anterior chamber
%   depth, refractive state, corneal diameter, and axial length."
%   Journal of Refractive Surgery 16.3 (2000): 336-340.
%
% To position the retina, we need to know the distance between the apex of
% the anterior chamber and the apex of the retina. I derive the value for
% this distance from the Atchison 2006 model eye.
retinaCenter = ...
    [(-(23.5800 - retinaRadiiEmetrope(1)*2) - retinaRadii(1)) 0 0];

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

