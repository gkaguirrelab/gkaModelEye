function [opticalSystemOut, p] = addContactLens(opticalSystemIn, lensRefractionDiopters, varargin)
% Add a contact lens to a passed optical system
%
% Syntax:
%  [opticalSystemOut, p] = addContactLens(opticalSystemIn, lensRefractionDiopters)
%
% Description:
%	This routine adds a meniscus (ophthalmologic) contact lens to an
%	optical system with the refractive power specified in the passed
%	variable. Note that a ray emerging from the eye encounters two concave
%	surfaces for this lens, so both surfaces will have a negative radius of
%	curvature for rayTraceCenteredSphericalSurfaces().
%
% Inputs:
%   opticalSystemIn       - An mx4 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           position of the ray. Each row contains the
%                           values:
%                               [center, radiusZ, radiusH, refractiveIndex]
%                           that define an elliptical lens.
%   lensRefractionDiopters - Scalar. Refractive power in units of
%                           diopters. A negative value specifies a lens
%                           that would be worn by someone with myopia to
%                           correct their vision.
%
% Optional key/value pairs:
%  'lensRefractiveIndex'  - Scalar. Refractive index of the lens material.
%                           The routine returnRefractiveIndex() provides
%                           the indices for several spectacle materials
%                           under visible (vis) and near infra-red (nir)
%                           imaging domains.
%
% Outputs:
%   opticalSystemOut      - An (m+1)x4 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           contact lens.
%   p                     - The parameters returned by the input parser.
%
% Examples:
%{
    %% Replicate calculation of contact lens curvature
    % WJ Benajamin provudes an example calculation for a contact lens
    % that provides -10D power. We confirm here that our routine provides
    % the same solution.
    %   Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact 
    %   lens practice. Lippincott Williams & Wilkins, 2005. Chapter 7A, 
    %   "Optical phenomena of contact lenses", WJ Benjamin. p130
    opticalSystemIn = [nan nan nan nan 1.3760; -7.8  -7.8  -7.8  -7.8  1.0];
    % The curvature of the front surface of the contact lens should be
    % -9.56. We obtain a slightly lower value (9.5468) as we do not model
    % the effect of the pre-lens tear film.
    opticalSystemOut = addContactLens(opticalSystemIn, -10, 'lensRefractiveIndex', 1.43 );
    assert(abs(opticalSystemOut(end,2) - -9.56)<0.1);
%}
%{
    %% Contact lens added to correct myopia
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
    rayTraceEllipsoids([sceneGeometry.eye.pupil.center(1),0],deg2rad(-15),sceneGeometry.refraction.opticalSystem,true);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('lensRefractionDiopters',@isnumeric);

% Optional
p.addParameter('lensRefractiveIndex',returnRefractiveIndex( 'hydrogel', 'NIR' ),@isnumeric);
p.addParameter('minimumLensThickness',0.05,@isnumeric);

% parse
p.parse(opticalSystemIn, lensRefractionDiopters, varargin{:})


% Distribute the parameters into variables
lensRefractiveIndex = p.Results.lensRefractiveIndex;

% Copy the optical system from input to output
opticalSystemOut = opticalSystemIn;

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. Because the contact lens contacts the
% corneal surface, this index of refraction will be replaced with the index
% of refraction of the contact lens. This is now the back surface of the
% contact lens. We store the index of refraction of the ambient medium
% (which will typically be air and thus 1.0) to apply to the final exit
% ray.
priorRefractiveIndex = opticalSystemIn(end-1,end);
mediumRefractiveIndex = opticalSystemIn(end,end);
opticalSystemOut(end,end) = lensRefractiveIndex;

% Calculate the diopters of the corneal surface without a contact lens. We
% consider only the radius of curvature at the apex along the optical axis.
% Our goal is to create a front surface of the contact lens that produces a
% refractive correction equal to:
%   cornealSurfaceDiopters + lensRefractionDiopters
t =0;
cornealSurfaceCurvature = -((opticalSystemIn(end,2)^2*sin(t)^2 + opticalSystemIn(end,3)^2*cos(t)^2)^(3/2))/(opticalSystemIn(end,2)*opticalSystemIn(end,3));
cornealSurfaceDiopters = (mediumRefractiveIndex-priorRefractiveIndex)/(cornealSurfaceCurvature/1000);

% Calculate the refractive power of the back surface of the lens.
backCurvature = cornealSurfaceCurvature;
backDiopters = (lensRefractiveIndex-priorRefractiveIndex)/(backCurvature/1000);

% We calculate here thickness and thus center of the front surface.
if lensRefractionDiopters > 0
    % This is a plus lens for the correction of hyperopia.
    
    % How many diopters of correction do we need from the front surface?
    frontDiopters = cornealSurfaceDiopters + lensRefractionDiopters - backDiopters;
    
    % Calculate the radius of curvature of the front surface. This is the
    % lens-maker's equation, re-arranged to solve for the initial radius of
    % curvature.
    %{
    syms r1 r2 n m t d
    lensMakersEqn = (n-m)*(1/r1 - 1/r2 + ((n-m)*t)/(n*r1*r2)) == d;
    frontCurvatureEqn = isolate(lensMakersEqn,r2)*1000;
    % where r1 = back curvature in meters; r2 = front curvature; n = index
    % of refraction of lens material; m = index of refraction of medium
    % (usually the air); t = thickness of the lens at the optical axis in
    % meters; d = diopters of the lens
    %}
    % First perform the calculation using the thin-lens approximation
    frontCurvatureApprox = (((0/1000*(mediumRefractiveIndex-lensRefractiveIndex))/ ...
        (lensRefractiveIndex*backCurvature/1000)+1) / ...
        (frontDiopters/(mediumRefractiveIndex-lensRefractiveIndex) + 1/backCurvature/1000))*1000;

    % Use the approximation to define a thickness at the optical axis
    thickness = max([frontCurvatureApprox-backCurvature p.Results.minimumLensThickness]);

    % Now calculate using the thickness
    frontCurvature = (((thickness/1000*(mediumRefractiveIndex-lensRefractiveIndex))/ ...
        (lensRefractiveIndex*backCurvature/1000)+1) / ...
        (frontDiopters/(mediumRefractiveIndex-lensRefractiveIndex) + 1/backCurvature/1000))*1000;

    % Calculate the location of the center of curvature for the front lens.
    frontCenter = frontCurvature + (frontCurvature-backCurvature);
    
else
    % This is a minus lens for the correction of myopia.
    % It will be thinnest at the center of the lens on the optical axis.
    % The parameter minimumLensThickness defines this minimum.
    
    % How many diopters of correction do we need from the lens?
    frontDiopters = cornealSurfaceDiopters + lensRefractionDiopters - backDiopters;
    
    % Calculate the radius of curvature of the front surface. This is the
    % lens-maker's equation, re-arranged to solve for the initial radius of
    % curvature.
    %{
    syms r1 r2 n m t d
    lensMakersEqn = (n-m)*(1/r1 - 1/r2 + ((n-m)*t)/(n*r1*r2)) == d;
    frontCurvatureEqn = isolate(lensMakersEqn,r2)*1000;
    % where r1 = back curvature in meters; r2 = front curvature; n = index
    % of refraction of lens material; m = index of refraction of medium
    % (usually the air); t = thickness of the lens at the optical axis in
    % meters; d = diopters of the lens
    %}
    frontCurvature = (((p.Results.minimumLensThickness/1000*(mediumRefractiveIndex-lensRefractiveIndex))/ ...
        (lensRefractiveIndex*backCurvature/1000)+1) / ...
        (frontDiopters/(mediumRefractiveIndex-lensRefractiveIndex) + 1/backCurvature/1000))*1000;
    
    % Calculate the location of the center of curvature for the front lens.
    frontCenter = frontCurvature + p.Results.minimumLensThickness;

end

% Add the contact lens to the optical system
opticalSystemOut(end+1,:)=[frontCenter frontCurvature frontCurvature frontCurvature mediumRefractiveIndex];

end % function - addContactLens