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
%   opticalSystemIn       - An mx3 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           position of the ray. Each row contains the
%                           values [center, radius, refractiveIndex] that
%                           define a spherical lens.
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
%   opticalSystemOut      - An (m+2)x3 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           contact lens
%   p                     - The parameters returned by the input parser.
%
% Examples:
%{
    %% Example - Replicate calculation of contact lens curvature
    % WJ Benajamin provudes an example calculation for a contact lens
    % that provides -10D power. We confirm here that our routine provides
    % the same solution.
    %   Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact 
    %   lens practice. Lippincott Williams & Wilkins, 2005. Chapter 7A, 
    %   "Optical phenomena of contact lenses", WJ Benjamin. p130
    eye = modelEyeParameters();
    eye.corneaFrontSurfaceR = 7.8;
    cornealThickness = eye.corneaBackSurfaceCenter(1) - eye.corneaFrontSurfaceCenter(1);
    opticalSystem = [nan, nan, eye.aqueousRefractiveIndex; ...
        -eye.corneaBackSurfaceR-cornealThickness, -eye.corneaBackSurfaceR, eye.corneaRefractiveIndex; ...
        -eye.corneaFrontSurfaceR, -eye.corneaFrontSurfaceR, 1.0];
    % Add a -10 diopter lens
    opticalSystem=addContactLens(opticalSystem, -10, 'lensRefractiveIndex', 1.43, 'minimumLensThickness', 0.1);
    % The curvature of the front surface of the contact lens should be
    % 9.56 mm
    opticalSystem(end,2)
%}
%{
    %% Example - Test the output for zero diopters
    % If zero diopters are requested for a contact lens with an index of
    % refraction equal to the corneal index, then the curvature of the
    % front surface of the lens should be equivalent to the corneal front
    % surface curvature.
    eye = modelEyeParameters();
    eye.corneaRefractiveIndex = 1.376;
    cornealThickness = eye.corneaBackSurfaceCenter(1) - eye.corneaFrontSurfaceCenter(1);
    opticalSystem = [nan, nan, eye.aqueousRefractiveIndex; ...
        -eye.corneaBackSurfaceR-cornealThickness, -eye.corneaBackSurfaceR, eye.corneaRefractiveIndex; ...
        -eye.corneaFrontSurfaceR, -eye.corneaFrontSurfaceR, 1.0];
    opticalSystem=addContactLens(opticalSystem, 0, 'lensRefractiveIndex', eye.corneaRefractiveIndex, 'minimumLensThickness',0);
    % Is the lens surface the same curvature as the corneal surface?
    round(opticalSystem(end,1),4) == round(opticalSystem(end-1,1),4)
%}
%{
    %% Example - Ray trace through cornea and contact lens
    % Obtain the eye parameters from the modelEyeParameters() function
    eye = modelEyeParameters('sphericalAmetropia',2);
    % Define an optical system
    cornealThickness = eye.corneaBackSurfaceCenter(1) - eye.corneaFrontSurfaceCenter(1);
    opticalSystem = [nan, nan, eye.aqueousRefractiveIndex; ...
        -eye.corneaBackSurfaceR-cornealThickness, -eye.corneaBackSurfaceR, eye.corneaRefractiveIndex; ...
        -eye.corneaFrontSurfaceR, -eye.corneaFrontSurfaceR, 1.0];
    % Add a plus lens for the correction of hyperopia
    opticalSystem=addContactLens(opticalSystem, 2);
    % Define FigureFlag as a structure so we can provide plot limits. Also,
    % turn off the text labels to reduce clutter
    clear figureFlag
    figureFlag.textLabels = false;
    figureFlag.zLim = [-10 5];
    figureFlag.hLim = [-5 5];
    % Define a ray originating from the border of the pupil
    pupilRadius = 2;
    clear coords
    clear theta
    theta = deg2rad(-30);
    coords = [eye.pupilCenter(1) pupilRadius];
    % Perform the ray tracing
    outputRay = rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem, figureFlag)
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
% contact lens.
% We store the index of refraction of the ambient medium (which will
% typically be air and thus 1.0) to apply to the final exit ray.
priorRefractiveIndex = opticalSystemIn(end-1,3);
mediumRefractiveIndex = opticalSystemIn(end,3);
opticalSystemOut(end,3) = lensRefractiveIndex;

% Calculate the diopters of the corneal surface without a contact lens; our
% goal is to create a front surface of the contact lens that produces a
% refractive correction equal to:
%   cornealSurfaceDiopters + lensRefractionDiopters
cornealSurfaceDiopters = (mediumRefractiveIndex-priorRefractiveIndex)/(opticalSystemIn(end,2)/1000);

% Calculate the refractive power of the back surface of the lens.
backCurvature = opticalSystemIn(end,2);
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
    
    % Store the lens front surface in the optical system
    opticalSystemOut(end+1,:)=[frontCenter frontCurvature mediumRefractiveIndex];
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

    % Add the surfaces to the optical system
    opticalSystemOut(end+1,:)=[frontCenter frontCurvature mediumRefractiveIndex];
end

end % function - addContactLens