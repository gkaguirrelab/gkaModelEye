function vitreousChamber = vitreousChamber(eye)

% Atchison 2006 provides radii of curvature and asphericities for a
% biconic model of the vitreous chamber, with these values varying
% by spherical ametropia. Parameters for the the decentration and
% tilt of the vitreous chamber are also provided:
%
%	Atchison, David A., et al. "Shape of the retinal surface in
%   emmetropia and myopia." Investigative ophthalmology & visual
%   science 46.8 (2005): 2698-2707.
%
% I model the vitreous chamber as a centered ellipsoid. I convert
% the 4 parameeters of the Atchison biconic model to a 3 radii of
% an ellipsoid by numeric approximation. To match Atchison's axial
% length formula (Eq 19), I had to inflate the effect of spherical
% ametropia upon the asphericity coefficients very sligtly.
% Atchison gives the values:
%
%   Qx = 0.27+0.026*SR
%   Qy = 0.25+0.017*SR
%
% and I increased the effects of SR to be 0.0272 and 0.0182 upon Qx
% and Qy, respectively. I suspect this adjustment is the result of
% a small, systematic underestimation of the ellipsoid radii by my
% numeric approximation.
%
%{
    % Numeric approximation of Atchison 2006 biconic model of
    % vitreous chamber with ellipsoid radii
    radii = [];
    for SR = -2:2
    Cx = 1/(12.91+0.094*SR);
    Cy = 1/(12.72-0.004*SR);
    Qx = 0.27+0.0272*SR;
    Qy = 0.25+0.0182*SR;
    biconicZ = @(x,y) (Cx.*x.^2 + Cy.*y.^2)./(1+sqrt( 1-(1+Qx).* Cx.^2.*x.^2 - (1+Qy).*Cy.^2.*y.^2));
    myObj = @(p) -biconicZ(p,0);
    [radiusX] = fminsearch(myObj,10);
    myObj = @(p) -biconicZ(0,p);
    [radiusY] = fminsearch(myObj,10);
    radiusZ = max([biconicZ(radiusX,0) biconicZ(0,radiusY)]);
    radii = [radii; [radiusZ, radiusX, radiusY]];
    end
    slopes = mean(diff(radii));
    fprintf('axial radius = %4.4f %4.4f * SR\n',radii(3,1),slopes(1));
    fprintf('horizontal radius = %4.4f %4.4f * SR \n',radii(3,2),slopes(2));
    fprintf('vertical radius = %4.4f %4.4f * SR \n',radii(3,3),slopes(3));
%}
vitreousChamberRadiiEmetrope = [10.1760 11.4558 11.3771];
vitreousChamberRadiiAmetropiaSlope = [-0.1495 -0.0393 -0.0864];
vitreousChamber.radii = ...
    vitreousChamberRadiiEmetrope + vitreousChamberRadiiAmetropiaSlope.* eye.meta.sphericalAmetropia;

% Our model holds the depth of the anterior chamber constant.
% Atchison found that anterior chamber depth does not vary with
% spherical ametropia, although this is not a consistent finding:
%
%   Hosny, Mohamed, et al. "Relationship between anterior chamber
%   depth, refractive state, corneal diameter, and axial length."
%   Journal of Refractive Surgery 16.3 (2000): 336-340.
%
% To position the vitreous chamber, we need to know the distance
% between the apex of the anterior chamber and the apex of the
% vitreous chamber. I derive the value for this distance from the
% Atchison 2006 model eye.
vitreousChamberApexDepth = 23.5800 - vitreousChamberRadiiEmetrope(1)*2;

% Set the depth of the center of the vitreous chamber
vitreousChamber.center = ...
    [(-vitreousChamberApexDepth - vitreousChamber.radii(1)) 0 0];

S = quadric.scale(quadric.unitSphere,vitreousChamber.radii);
S = quadric.translate(S,vitreousChamber.center);
vitreousChamber.S = quadric.matrixToVec(S);
vitreousChamber.side = -1;
vitreousChamber.boundingBox = [-25 -5.4 -25 25 -25 25];
vitreousChamber.mustIntersect = 1;
vitreousChamber.label = {'vitreousChamber'};
vitreousChamber.plot.color = {[.7,.5,.7]};

end
