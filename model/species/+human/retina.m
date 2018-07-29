function retina = retina(eye)

% Atchison 2006 provides the semi-radii of a triaxial ellipsoid model of
% the retinal surface, with these values varying by spherical ametropia:
%
%	Atchison, David A., et al. "Shape of the retinal surface in
%   emmetropia and myopia." Investigative ophthalmology & visual
%   science 46.8 (2005): 2698-2707.
%
% The Atchison model is arranged around the visual axis, with the vitreous
% chamber tilted and shifted. In the current model, the retinal ellipsoid
% is aligned with and centered on the optical axis.

retinaRadiiEmetrope = [10.148 11.455 11.365];
retinaRadiiAmetropiaSlope = [-0.163 -0.043 -0.090];
retinaRadii = ...
    retinaRadiiEmetrope + retinaRadiiAmetropiaSlope.* eye.meta.sphericalAmetropia;

% Our model holds the depth of the anterior chamber constant.
% Atchison found that anterior chamber depth does not vary with
% spherical ametropia, although this is not a consistent finding:
%
%   Hosny, Mohamed, et al. "Relationship between anterior chamber
%   depth, refractive state, corneal diameter, and axial length."
%   Journal of Refractive Surgery 16.3 (2000): 336-340.
%
% To position the retina, we need to know the distance
% between the apex of the anterior chamber and the apex of the
% retina. I derive the value for this distance from the
% Atchison 2006 model eye.
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

