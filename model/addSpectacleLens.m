function [opticalSystemOut, p] = addSpectacleLens(opticalSystemIn, lensRefractionDiopters, varargin)
% Add a spectacle lens to a passed optical system
%
% Syntax:
%  [opticalSystemOut, p] = addSpectacleLens(opticalSystemIn, lensRefractionDiopters)
%
% Description:
%	This routine adds a meniscus (ophthalmologic) spectacle lens to an
%	optical system with the refractive power specified in the passed
%	variable. Note that a ray emerging from the eye encounters two concave
%	surfaces for this lens, so both surfaces will have a negative radius of
%	curvature for rayTraceEllipsoids().
%
% Inputs:
%   opticalSystemIn       - An mx19 matrix, where m is set by the key value
%                           opticalSystemNumRows. Each row contains the 
%                           values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
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
%  'lensVertexDistance'   - Scalar. Distance (in mm) between the corneal
%                           apexa and the back surface of the lens. Typical
%                           values are 12-14 mm.
%  'nearPlanoCurvature'   - Scalar. This defines the curvature of the near-
%                           plano face of the lens.
%
% Outputs:
%   opticalSystemOut      - An (m+2)x19 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           spectacle lens.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('lensRefractionDiopters',@isnumeric);

% Optional
p.addParameter('lensRefractiveIndex',returnRefractiveIndex( 'polycarbonate', 'NIR' ),@isnumeric);
p.addParameter('lensVertexDistance',12,@isnumeric);
p.addParameter('nearPlanoCurvature',-16000,@isnumeric);
p.addParameter('minimumLensThickness',1,@isnumeric);

% parse
p.parse(opticalSystemIn, lensRefractionDiopters, varargin{:})

% Distribute the parameters into variables
lensVertexDistance = p.Results.lensVertexDistance;
lensRefractiveIndex = p.Results.lensRefractiveIndex;
nearPlanoCurvature = p.Results.nearPlanoCurvature;

% Obtain the quadric for the corneal surface
corneaLine = opticalSystemIn(end,:);

% Copy the optical system from input to output
opticalSystemOut = opticalSystemIn;

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.
mediumRefractiveIndex = opticalSystemIn(end,end);

% The lens equations do not perform properly for corrections of less that
% 0.25 diopters, and we don't bother trying to model so small a correction.
% In such a case, return the opticalSystem unaltered.
if abs(lensRefractionDiopters) < 0.25
    return
end

% Define a "meniscus" lens with two surfaces.
if lensRefractionDiopters > 0
    % This is a plus lens for the correction of hyperopia. It has a
    % relatively flat back surface and a more curved front surface.
    % We first add the near-plano back surface to the optical system
    backCurvature = nearPlanoCurvature;
    backCenter = lensVertexDistance+backCurvature;
    backDiopters = (mediumRefractiveIndex-lensRefractiveIndex)/(backCurvature/1000);
 
    % How many diopters of correction do we need from the front surface?
    frontDiopters = lensRefractionDiopters - backDiopters;
 
    % Calculate the radius of curvature of the front surface using the thin
    % lens approximation
    frontCurvature = ((mediumRefractiveIndex-lensRefractiveIndex)/frontDiopters)*1000;
    
    % Calculate the location of the center of curvature for the front lens.
    % To do so, we introduce lens thickness as a symbolic variable.
    syms thickness
    assume(thickness>0);
    frontCenter = lensVertexDistance + frontCurvature + thickness;
        
    % Define the equations for the height of each surface
    syms x
    assume(x>0);
    backPerimHeight = sqrt(backCurvature^2 - (x-backCenter)^2);
    frontPerimHeight = sqrt(frontCurvature^2 - (x-frontCenter)^2);
    
    % The lens will be thickest along the optical axis of the eye. We wish
    % the lens to have a non-negative thickness within the range of rays
    % we wish to model as coming from the eye. We consider a line with a
    % slope of 2 that originates at the corneal apex and find the height at
    % which it intersects the back surface.
    eqn= backPerimHeight == x*2;
    intersectHeight = eval(solve(eqn));

    % Now solve for the thickness
    eqn = subs(frontPerimHeight,x,intersectHeight) == subs(backPerimHeight,x,intersectHeight);
    thickness = min(eval(solve(eqn)));

    % Clear the symbolic variables
    clear x

    frontCenter = double(subs(frontCenter));
        
else
    % This is a minus lens for the correction of myopia. It has a
    % relatively flat front surface and a more curved back surface. It will
    % be thinnest at the center of the lens on the optical axis. The
    % parameter minimumLensThickness defines this minimum.

    % We first determine the properties of the the near-plano front surface
    frontCurvature = nearPlanoCurvature;
    frontDiopters = (mediumRefractiveIndex-lensRefractiveIndex)/(frontCurvature/1000);
 
    % How many diopters of correction do we need from the back surface?
    backDiopters = frontDiopters - lensRefractionDiopters;
 
    % Calculate the radius of curvature of the back surface using the thin
    % lens approximation
    backCurvature = ((mediumRefractiveIndex-lensRefractiveIndex)/backDiopters)*1000;
    
    % Calculate the locations of the center of curvature for the back and
    % front lens.
    backCenter = lensVertexDistance+backCurvature;
    frontCenter = lensVertexDistance+frontCurvature+p.Results.minimumLensThickness;
    
end % positive or negative lens

% Define a bounding box
boundingBoxLens = [0 frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

% Add the back spectacle surface to the optical system.
SlensBack = quadric.scale(quadric.unitSphere,[backCurvature backCurvature backCurvature]);
SlensBack = quadric.translate(SlensBack,[backCenter 0 0]);
lensLine = corneaLine;
lensLine(1:10) = quadric.matrixToVec(SlensBack);
lensLine(12:17) = boundingBoxLens;
lensLine(end) = p.Results.lensRefractiveIndex;
opticalSystemOut = [opticalSystemOut; lensLine];

% Add the back spectacle surface to the optical system.
SlensFront = quadric.scale(quadric.unitSphere,[frontCurvature frontCurvature frontCurvature]);
SlensFront = quadric.translate(SlensFront,[frontCenter 0 0]);
lensLine = corneaLine;
lensLine(1:10) = quadric.matrixToVec(SlensFront);
lensLine(12:17) = boundingBoxLens;
lensLine(end) = mediumRefractiveIndex;
opticalSystemOut = [opticalSystemOut; lensLine];

end % function - addSpectacleLens