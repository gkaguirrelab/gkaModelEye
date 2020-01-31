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
% Examples:
%{
    % Navarro 2006 examines the orientation of the corneal ellipsoid with
    % respect to the keratometric axis, which is the line that connects the 
    % fixation point of a keratometer instrument with the corneal center of
    % curvature. Here I obtain the distance in mm between the modeled apex 
    % of the cornea, and the "vertex normal" point, which is the point of
    % intersection of the keratometric axis with the front surface of the
    % cornea. The calculation is made for an emmetropic eye and assumes 
    % that the fixation point of the keratometric instrument is at 500 mm  
    % from the corneal surface.

    % Create a sceneGeometry. Put the fixation target at 500 mm, set the
    % stop radius to 0.86 mm, which corresponds to a 2 mm diameter pupil.

	fixTargetDistance = 500;
	stopRadius = 0.8693;
    sceneGeometry = createSceneGeometry(...
        'sphericalAmetropia',0,...
        'accommodationDiopters',1000/fixTargetDistance,...
        'spectralDomain','vis',...
        'calcLandmarkFovea',true);

    % Obtain the fixation angles and fixation target location
	[~,~,~, fixTargetWorldCoords] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance);

	% Find the point of intersection of the keratometric ray with the
	% cornea
	fixTargetEyeCoords = fixTargetWorldCoords([3 1 2]);
	cc = sceneGeometry.eye.cornea.front.center;
	keratometricAxisRay = quadric.normalizeRay([fixTargetEyeCoords';cc-fixTargetEyeCoords']');
    opticalSystem = sceneGeometry.refraction.cameraToRetina.opticalSystem;
    [~,rayPath] = rayTraceQuadrics(keratometricAxisRay, opticalSystem);
    horizontalDistance = rayPath(2,2);
%}


%% Corneal ellipsoid rotation
% The corneal ellipsoid is modeled as aligned with the optical axis. The
% code below retains the ability to model a rotation if desired. The
% cornealRotation vector specifies the rotation (in degrees) about each of
% the axes.
cornealRotation = [0 0 0];

%% Front corneal surface
if isempty(eye.meta.measuredCornealCurvature)
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
    % We calculate the change in parameters of the Navarro model that would
    % be expected given the Atchison effect for ametropia.
    %{
        R = @(D) 7.77 + 0.022 .* D;
        Q = -0.15;
        a = @(D) R(D) ./ (Q+1);
        b = @(D) R(D) .* sqrt(1./(Q+1));
        radiiAtchFront = @(D) [a(D) b(D) b(D)];
        % Show that the ametropia correction scales all radii equally
        radiiAtchFront(0)./radiiAtchFront(1)
        % Calculate the proportion change in radius
        radiusScalerPerD = 1-a(1)/a(0);
        radiiNavFront = [14.26   10.43   10.27];
        radiiNavFrontCorrected = @(D) radiiNavFront.* (D.*radiusScalerPerD+1);
        % Report the ratio of the Atchison and Navarro axial radii for the
        % front surface of the cornea; we use this below.
        atchNavScaler = a(0) ./ radiiNavFront(1)
    %}
    radii = [14.26   10.43   10.27] .* ...
        ((eye.meta.sphericalAmetropia .* -0.0028)+1);
    S = quadric.scale(quadric.unitSphere,radii);    
else
    % If a measured value is provided, use it here to calculate the
    % parameters of the ellipsoidal surface. We set the axial length of the
    % corneal ellipsoid equal to the value provided by Navarro 2006.
    radii(1) = 14.26;
    % Formula to convert diopters to radius of curvature in mm
    RoC = @(D) 1000.*(1.3375-1)./D;
    % The horizontal and vertical radii are derived from the passed values
    radii(2:3) = sqrt(radii(1).*RoC(eye.meta.measuredCornealCurvature(1:2)));
    % Create the quadric
    S = quadric.scale(quadric.unitSphere,radii);
    % Here is where I would apply a torsional rotation to the ellipsoid.
    % However, the quadric rotation is not working properly so this is not
    % implemented currently.
    if length(eye.meta.measuredCornealCurvature)==3
%         cornealRotation = [eye.meta.measuredCornealCurvature(3) 0 0];
    end
end

% Saving this for use in constructing the tear film below
S_front = S;

% Rotate the quadric surface
switch eye.meta.eyeLaterality
    case 'Right'
        S = quadric.rotate(S,cornealRotation);
    case 'Left'
        S = quadric.rotate(S,cornealRotation.*[1 1 -1]);
    otherwise
        error('eye laterality not defined')
end

% We set the center of the cornea front surface ellipsoid so that the axial
% apex (prior to rotation) is at position [0, 0, 0]
S = quadric.translate(S,[-radii(1) 0 0]);

% Store these values
cornea.front.S = quadric.matrixToVec(S);
cornea.front.side = 1;
cornea.front.boundingBox=[-4 0 -8 8 -8 8];
cornea.front.center=[-radii(1) 0 0];


%% Tear film
% The tear film is the front corneal surface, translated forward, and with
% a different refractive index. The thickness is taken from:
%   Werkmeister, René M., et al. "Measurement of tear film thickness using
%   ultrahigh-resolution optical coherence tomography." Investigative
%   ophthalmology & visual science 54.8 (2013): 5578-5583.
tearFilmThickness = 0.005;

% Rotate the quadric surface
S = S_front;
switch eye.meta.eyeLaterality
    case 'Right'
        S = quadric.rotate(S,cornealRotation);
    case 'Left'
        S = quadric.rotate(S,cornealRotation.*[1 1 -1]);
    otherwise
        error('eye laterality not defined')
end

% Translate and store
S = quadric.translate(S,[-radii(1)+tearFilmThickness 0 0]);
cornea.tears.S = quadric.matrixToVec(S);
cornea.tears.side = 1;
cornea.tears.boundingBox=[-4+tearFilmThickness tearFilmThickness -8 8 -8 8];


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
radii = [ 13.7716    9.3027    9.3027];
S = quadric.scale(quadric.unitSphere,radii);

% Rotate the quadric surface
switch eye.meta.eyeLaterality
    case 'Right'
        S = quadric.rotate(S,cornealRotation+[ 180 180 180 ]);
    case 'Left'
        S = quadric.rotate(S,cornealRotation.*[1 1 -1]-[ 180 180 180 ]);
    otherwise
        error('eye laterality not defined')
end

% The center of the back cornea ellipsoid is positioned so that there is
% 0.55 mm of corneal thickness between the front and back surface of the
% cornea at the apex, following Atchison 2006.
cornealThickness = 0.55;
S = quadric.translate(S,[-cornealThickness-radii(1) 0 0]);

% Store these values
cornea.back.S = quadric.matrixToVec(S);
cornea.back.side = 1;
cornea.back.boundingBox=[-4 0 -8 8 -8 8];

% Assemble the combined corneal surfaces
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