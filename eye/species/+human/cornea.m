function cornea = cornea( eye )
% Returns the cornea sub-field of an eye model structure
%
% Syntax:
%  cornea = human.cornea( eye )
%
% Description:
%   The corneal front surface is taken from Table 1 of Navarro 2006:
%
%       Navarro, Rafael, Luis González, and José L. Hernández. "Optics of
%       the average normal cornea from general and canonical
%       representations of its surface topography." JOSA A 23.2 (2006):
%       219-232.
%
%   Their dimensions [a,b,c] correspond to [p2, p3, p1] in the current
%   model. The cornea ellipsoid is modeled as being aligned with the
%   optical axis.
%
%   The radius of curvature at the vertex of the cornea was found by
%   Atchison to vary as a function of spherical ametropia (Table 1):
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision
%       research 46.14 (2006): 2236-2250.
%
%   The Navarro parameters are adjusted here to reflect this variation by
%   ametropia.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   cornea                - Structure.
%


%% Cornea axial radius
% The front surface of the cornea is well modeled as a tri-axial ellipsoid.
% Measurements of the cornea are frequently specified as k1 and k2 values
% (which are in units of diopters) or in terms of the radius of curvature
% (R) at the vertex and asphericity (Q). In both of these frameworks, it is
% also necessary to specify the radius of the corneal ellipsoid in the
% axial direction.
%
% Navarro 2006 provides a value for the axial radius of 14.26 mm (for a
% myopic population). Atchison 2006 observes that the radius of curvature
% (but not the asphericity) of the cornea front surface varies with
% spherical ametropia. This variation amounts to a change in the axial
% radius (without a change in refractive power). If not otherwise
% specified, the cornea axial radius is set to the Navarro value, adjusted
% for the ametropia effect reported by Atchison.

if isempty(eye.meta.corneaAxialRadius)
                
    % Atchison provides parameters for a radially symmetric ellipsoid in
    % terms of the radius of curvature (R) at the vertex and its
    % asphericity (Q). R varies with spherical ametropia (D):
    %
    %   R = 7.77 + 0.022 * D
    %   Q = -0.15
    %
    % Because the asphericity of the cornea did not change, the change in R
    % corresponds to an overall scaling of the ellipsoid in all dimensions.
    % We adjust the Navarro values to account for this effect. R and Q are
    % related to the radii of an ellipse along the primary and secondy axes
    % (a, b) by:
    %
    %   R = b^2/a
    %	Q = (b^2 / a^2) - 1
    %
    % when Q < 0. Therefore, given R and Q, we can obtain a and b, which
    % correspond to the radii of the ellipsoid model, with a corresponding
    % to the axial dimension, and b to the horizontal and vertical
    % dimensions. Checking my algebra here:
    %{
        syms a b R Q
        eqn1 = R == b^2/a;
        eqn2 = Q == (b^2 / a^2) - 1;
        solution = solve([eqn1, eqn2]);
        solution.a
        solution.b
    %}
    % We calculate the multiplier effect on corneal axial radius implied by
    % the Atchison effect
    %{
        R = @(D) 7.77 + 0.022 .* D;
        Q = -0.15;
        a = @(D) R(D) ./ (Q+1);
        b = @(D) R(D) .* sqrt(1./(Q+1));
        % Calculate the proportion change in radius
        radiusScalerPerD = 1-a(1)/a(0);
    %}
    radiusScalerPerD = -0.002831402831403;
    
    % Navarro 2006 provides the mean cornea axial radius for a population
    % of people with an average spherical refractive error of -3.9 D. 
    navarroMyopicRadius = 14.26;
    navarroPopulationAmetropia = -3.9;
    
    % Using the Atchison result, we can back out the effect of myopia in
    % the Navarro population
    corneaAxialRadiusEmmetrope = navarroMyopicRadius / ...
        ((navarroPopulationAmetropia * radiusScalerPerD)+1);
    
    % Finally, we add back in the effect of ametropia in the eye to be
    % modeled
    corneaAxialRadius = corneaAxialRadiusEmmetrope * ...
        ((eye.meta.sphericalAmetropia * radiusScalerPerD)+1);
    
else
    
    % If the axial radius of the cornea was explicitly specified, use that
    % value.   
    corneaAxialRadius = eye.meta.corneaAxialRadius;

end

% Store the value
cornea.axialRadius = corneaAxialRadius;


%% Diopters to transverse radii
% The front of the cornea is a refractive surface, the optical power of
% which can be specified in units of keratometric diopters. Because it is
% aspheric, two values are used. The first value is always smaller, and
% thus describes the "flatter" surface of the cornea. The second value is
% always larger, and describes the curvature of the surface of the cornea
% oriented 90 degrees away from the flatest surface.
%
% These k values are returned by an ophthalmologic instrument that measures
% corneal curvature (expressed as radius of curvature at the corneal apex)
% and then converts that number into an equivalent optical power. The
% optical power reported is tweaked a bit so that the effect of the
% posterior surface of the cornea is included, and to force a relationship
% that Roc 7.5 mm == 45 diopters. To do so, a "keratometric" index of
% refraction of 1.3375 for the cornea is assumed, yielding this equation:
%
%       k (diopters) = (1.3375 - 1) / RoC (meters)
%
% Given the radius-of-curvature, we can obtain the corresponding semi-axes
% of the corneal ellipsoid if we know the length of the semi-axis of the
% cornea in the axial direction.
%
% These anonymous functions convert between keratometric power (in
% diopters) to the corresponding ellipsoid semi-axis length (in mm).
%
radiusFromPower = @(k) sqrt(corneaAxialRadius.*1000.*(1.3375-1)./k);
powerFromRadius = @(r) (corneaAxialRadius * 337.5) ./ r.^2;


%% Front surface power
% There are multiple studies of the distribution of normal corneal biometry
% values. These studies report effects of age, ametropia, and gender. If
% not otherwise specified, we adopt here the mean k values obtained from
% the normal population in the ConesToCortex human connectome project
if isempty(eye.meta.kvals)
    cornea.kvals(1:2) = [43.399, 44.33653846];
else
    cornea.kvals = eye.meta.kvals;
end
frontSurfaceRadii = ...
    [corneaAxialRadius, radiusFromPower(cornea.kvals(1)), radiusFromPower(cornea.kvals(2))];
cornea.frontSurfaceRadii = frontSurfaceRadii;

%% Cornea ellipsoid rotation
% The cornea can be rotated out of alignment with the optical axis of the
% eye, although the default value is to maintain the cornea centered on the
% longitudinal axis. If rotation were to be added, one option is to follow
% Navarro 2006 and rotate the apex of the corneal ellipsoid towards the
% visual axis of the eye, but not completely. About 2.5 degrees of "tip"
% rotation about the vertical axis towards the nose would be a reasonable
% choice.
%
% Note the order of listing of angles. In the kvals vector, the rotations
% are listed in the order:
%
%   [ torsion, tip (rotation about the horizontal axis, tilt (rotation about vertical) ]
%
% But in the vector that applies the rotation, we must give these angles in
% the order:
%
%   [torsion, tilt (rotation about vertical), tip (rotation about the horizontal axis) ]
%
cornea.kvals(3:5) = [0 0 0];
if ~isempty(eye.meta.kvals)
    nkvals = length(eye.meta.kvals);
    if nkvals > 2
        cornea.kvals(3:nkvals) = eye.meta.kvals(3:nkvals);
    end
end

% Create the corneal rotation vector. Note the different angle order from
% the kvals vector
corneaRotation = [cornea.kvals(3) cornea.kvals(5) cornea.kvals(4)];

% Account for the effects of eye laterality upon the angles
switch eye.meta.eyeLaterality
    case 'Right'
        % nothing to do
    case 'Left'
        corneaRotation = corneaRotation.*[1 1 -1];
    otherwise
        error('eye laterality not defined')
end

% Store the rotation vector in the cornea structure
cornea.rotation = corneaRotation;


%% Front corneal surface

% Create the quadric
S = quadric.scale(quadric.unitSphere,frontSurfaceRadii);

% Rotate the quadric surface
S = quadric.rotate(S,cornea.rotation);

% We set the center of the cornea front surface ellipsoid so that the axial
% apex (prior to rotation) is at position [0, 0, 0]
S = quadric.translate(S,[-frontSurfaceRadii(1) 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store these values
cornea.front.S = quadric.matrixToVec(S);
cornea.front.side = 1;
cornea.front.boundingBox=[-4 X(1) -9 9 -9 9];
cornea.front.center=[-frontSurfaceRadii(1) 0 0];


%% Tear film
% The tear film is the front corneal surface, translated forward, and with
% a different refractive index. The thickness is taken from:
%
%   Werkmeister, René M., et al. "Measurement of tear film thickness using
%   ultrahigh-resolution optical coherence tomography." Investigative
%   ophthalmology & visual science 54.8 (2013): 5578-5583.
tearFilmThickness = 0.005;

% Create the quadric
S = quadric.scale(quadric.unitSphere,frontSurfaceRadii);

% Rotate the quadric surface
S = quadric.rotate(S,cornea.rotation);

% Translate
S = quadric.translate(S,[-frontSurfaceRadii(1)+tearFilmThickness 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store
cornea.tears.S = quadric.matrixToVec(S);
cornea.tears.side = 1;
cornea.tears.boundingBox=[-4+tearFilmThickness X(1) -8 8 -8 8];


%% Back corneal surface
% Atchison finds that the back surface of the cornea does not vary by
% ametropia. Navarro does not provide posterior cornea parameters.
% Therefore, we scale the parameters provided by Atchison to relate to the
% axial corneal radius specified by Navarro:
%{
    R = 6.4;
    Q = -0.275;
    a = R ./ (Q+1);
    b = R .* sqrt(1./(Q+1));
    % Taken from the prior block of code
    atchNavScaler = 0.6410;
    radiiAtchBack = [a b b];
    % Scale the overall back cornea ellipsoid to match Navarro
    radiiNavBack = radiiAtchBack./atchNavScaler;
    % Now scale the relative horizontal and vertical axes so that
    % the relationship between the horizontal (and vertical) radii
    % and the axial radius is of the same proportion to the front
    % surface in the Atchison model
    radiiAtchFront0D = radiiAtchFront(0);
    frontHorizToAxAtch = radiiAtchFront0D(2)/radiiAtchFront0D(1);
    backHorizToAxAtch = b / a;
    radiiNavFront0D = radiiNavFrontCorrected(0);
    frontHorizToAxNav = radiiNavFront0D(2)/radiiNavFront0D(1);
    backHorizToAxNav = radiiNavBack(2)/radiiNavBack(1);
    targetBackHorizToAxNav = backHorizToAxAtch / frontHorizToAxAtch * frontHorizToAxNav;
    radiiNavBackCorrected = [a a*targetBackHorizToAxNav a*targetBackHorizToAxNav]./atchNavScaler
%}
backSurfaceRadii = [ 13.7716    9.3027    9.3027];
S = quadric.scale(quadric.unitSphere,backSurfaceRadii);

% Rotate the quadric surface
S = quadric.rotate(S,cornea.rotation);

% The center of the back cornea ellipsoid is positioned so that there is
% 0.55 mm of corneal thickness between the front and back surface of the
% cornea at the apex, following Atchison 2006.
cornealThickness = 0.55;
S = quadric.translate(S,[-cornealThickness-backSurfaceRadii(1) 0 0]);

% Find the most anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store these values
cornea.back.S = quadric.matrixToVec(S);
cornea.back.side = 1;
cornea.back.boundingBox=[-4 X(1) -9 9 -9 9];


%% Assemble the combined corneal surfaces
cornea.S = [cornea.back.S; cornea.front.S; cornea.tears.S];
cornea.boundingBox = [cornea.back.boundingBox; cornea.front.boundingBox; cornea.tears.boundingBox];
cornea.side = [1; 1; 1];
cornea.mustIntersect = [1; 1;1 ];
cornea.index = [returnRefractiveIndex( 'cornea', eye.meta.spectralDomain); ...
    returnRefractiveIndex( 'tears', eye.meta.spectralDomain)];
cornea.label = {'cornea.back'; 'cornea.front'; 'cornea.tearfilm'};
cornea.plot.color = {[0.5 0.5 0.75]; [0.5 0.5 0.75]; [0.5 0.5 0.75]};

% Code here to calculate the Navarro 1985 corneal parameters that
% were used by Fedtke 2010 in her simulation. These may be used for
% comparison.
%{
    % cornea front
    R = 7.72;
    Q = -0.26;
    a = R ./ (Q+1);
    b = R .* sqrt(1./(Q+1));
    [a b b]
    % cornea back
    R = 6.5;
    Q = 0;
    a = R ./ (Q+1);
    b = R .* sqrt(1./(Q+1));
    [a b b]
    cornea.front.radii = [10.4324    8.9743    8.9743];
    cornea.back.radii = [6.5000    6.5000    6.5000];
%}

end