function [opticalSystemOut, p, spectacleMagnification] = addSpectacleLens(opticalSystemIn, lensRefractionDiopters, varargin)
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
%                           This may be passed as empty, in which case the
%                           returned optical system will originate in air.
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
%  'minimumLensThickness' - Scalar. The minimum lens in mm.
%
% Outputs:
%   opticalSystemOut      - An (m+2)x19 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           spectacle lens.
%
% Examples:
%{
    % Confirm that a spectacle lens has the specified optical power
    lensVertexDistance = 14;
    entrancePupilDepth = 3;
    lensDiopters = 5;
    opticalSystem = addSpectacleLens([],lensDiopters,'lensVertexDistance',lensVertexDistance,'entrancePupilDepth',entrancePupilDepth,'systemDirection','cameraToEye');
    % Plot this
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true);
    % Trace parallel rays from right (the world) to left (the eye)
    R1 = quadric.normalizeRay([50,-1;-3,0;0,0]);
    R2 = quadric.normalizeRay([50,-1;3,0;0,0]);
    [outputRay1,rayPath1] = rayTraceQuadrics(R1, opticalSystem);
    [outputRay2,rayPath2] = rayTraceQuadrics(R2, opticalSystem);    
    % Add these rays to the plot
    plotOpticalSystem('newFigure',false,'outputRay',outputRay1,'rayPath',rayPath1);
    plotOpticalSystem('newFigure',false,'outputRay',outputRay2,'rayPath',rayPath2);
    % Calculate the focal point from the output rays. Calculate the lens
    % and compare to the called for value
    focalPoint=quadric.distanceRays(outputRay1,outputRay2);
    calcDiopters = -1000 / (focalPoint(1) - (lensVertexDistance+entrancePupilDepth));
    assert(abs((calcDiopters - lensDiopters)/lensDiopters) < 0.025);
%}
%{
    % Confirm that the returned spectacle magnification matches the
    % example given here:
    % http://www.drdrbill.com/downloads/optics/ophth-optics/Spectacle_Magnification.pdf
    % This is a +20 lens, with a +12 base curve, and a vertex distance from
    % the entrance pupil of 25 mm (3 mm added within the function for the
    % cornea to pupil distance). The given solution is x1.3 magnification.
    [~, ~, spectacleMagnification]=addSpectacleLens([],20,'lensRefractiveIndex',1.5,'lensVertexDistance',22,'baseCurve',12);
    assert(abs(spectacleMagnification-1.3)<0.1);
%}
%{
    % Calculate magnification by ray-tracing
    lensDiopters = -5;
    opticalSystem=addSpectacleLens([],lensDiopters,'systemDirection','cameraToEye');
    % Trace two rays from right (the world) to left (the eye) at different
    % initial angles
    height = 3;
    distance = 500;
    R1 = quadric.normalizeRay([distance,-1;height,0;0,0]);
    R2 = quadric.normalizeRay([distance,-1;height,-0.01;0,0]);
    [outputRay1,rayPath1] = rayTraceQuadrics(R1, opticalSystem);
    [outputRay2,rayPath2] = rayTraceQuadrics(R2, opticalSystem);    
    % Calculate the point of intersection of these two object rays
    imagePoint=quadric.distanceRays(outputRay1,outputRay2);
    imageHeight = imagePoint(2);
    magnification = imageHeight / height;
%}



%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('lensRefractionDiopters',@isnumeric);

% Optional
p.addParameter('lensRefractiveIndex',returnRefractiveIndex( 'polycarbonate', 'NIR' ),@isnumeric);
p.addParameter('lensVertexDistance',12,@isnumeric);
p.addParameter('baseCurve',[],@isnumeric);
p.addParameter('entrancePupilDepth',3,@isnumeric);
p.addParameter('minimumLensThickness',0.8,@isnumeric);
p.addParameter('systemDirection','eyeToCamera',@ischar);

% parse
p.parse(opticalSystemIn, lensRefractionDiopters, varargin{:})

% Distribute the parameters into variables
lensVertexDistance = p.Results.lensVertexDistance;
lensRefractiveIndex = p.Results.lensRefractiveIndex;

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.
if isempty(opticalSystemIn)
    mediumRefractiveIndex = 1;
else
    mediumRefractiveIndex = opticalSystemIn(end,end);
end

% The lens equations do not perform properly for corrections of less that
% 0.25 diopters, and we don't bother trying to model so small a correction.
% In such a case, return the opticalSystem unaltered.
if abs(lensRefractionDiopters) < 0.25
    spectacleMagnification = 1;
    return
end

% Set the base curve (front surface refraction in diopters) using Vogel's
% Rule if not otherwise specified
if isempty(p.Results.baseCurve)
    if lensRefractionDiopters > 0
        frontDiopters = lensRefractionDiopters + 6;
    else
        frontDiopters = lensRefractionDiopters/2 + 6;
    end
else
    frontDiopters = p.Results.baseCurve;
end

% Obtain the curvature of the front surface using the thin lens
% approximation
frontCurvature = (mediumRefractiveIndex-lensRefractiveIndex)/frontDiopters*1000;

% Define an "ophthalmic"  or convex-concave lens with two surfaces.
if lensRefractionDiopters > 0
    % This is a plus lens for the correction of hyperopia. How many
    % diopters of correction do we need from the back surface?
    backDiopters = lensRefractionDiopters - frontDiopters;
    
    % Calculate the radius of curvature of the back surface using the thin
    % lens approximation
    backCurvature = ((mediumRefractiveIndex-lensRefractiveIndex)/backDiopters)*1000;
    
    % The center of the back lens
    backCenter = lensVertexDistance - abs(backCurvature);
    
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
    % flatter front surface and a more curved back surface. It will
    % be thinnest at the center of the lens on the optical axis. The
    % parameter minimumLensThickness defines this minimum.
    thickness = p.Results.minimumLensThickness;
    
    % How many diopters of correction do we need from the back surface?
    backDiopters = lensRefractionDiopters - frontDiopters;
    
    % Calculate the radius of curvature of the back surface using the thin
    % lens approximation. Need a minus sign to properly position the lens.
    backCurvature = -((mediumRefractiveIndex-lensRefractiveIndex)/backDiopters)*1000;
    
    % Calculate the locations of the center of curvature for the back and
    % front lens.
    backCenter = lensVertexDistance+backCurvature;
    frontCenter = lensVertexDistance+frontCurvature+thickness;
    
    % Define the equations for the height of the back surface
    syms x
    assume(x>0);
    backPerimHeight = sqrt(backCurvature^2 - (x-backCenter)^2);
    
    % Consider a line with a slope of 2 that originates at the corneal apex
    % and find the height at which it intersects the back surface.
    eqn= backPerimHeight == x*2;
    intersectHeight = eval(solve(eqn));
    
end % positive or negative lens


%% Calculate the spectacle magnification
% The magnification produced by a thick spectacle lens has two components,
% the shape factor and the power factor
sf = 1 / (1 - ( (thickness/1000)/lensRefractiveIndex)*frontDiopters);
pf = 1 / (1 - (((p.Results.lensVertexDistance+p.Results.entrancePupilDepth)/1000)*backDiopters));
spectacleMagnification = (sf * pf);


%% Assemble the optical system
% If an optical system was provided, add to it. Otherwise, create an
% optical system that begins in air.
if isempty(opticalSystemIn)
    opticalSystemOut = nan(1,19);
    opticalSystemOut(19) = 1;
else
    % Copy the optical system from input to output
    opticalSystemOut = opticalSystemIn;
end

% Define a bounding box
boundingBoxLens = [intersectHeight frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

switch p.Results.systemDirection
    case 'eyeToCamera'
        % Add the back spectacle surface to the optical system.
        SlensBack = quadric.scale(quadric.unitSphere,[backCurvature backCurvature backCurvature]);
        SlensBack = quadric.translate(SlensBack,[backCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensBack);
        lensLine(11) = 1; % rays intersect concave lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = p.Results.lensRefractiveIndex;
        opticalSystemOut = [opticalSystemOut; lensLine];
        
        % Add the front spectacle surface to the optical system.
        SlensFront = quadric.scale(quadric.unitSphere,[frontCurvature frontCurvature frontCurvature]);
        SlensFront = quadric.translate(SlensFront,[frontCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensFront);
        lensLine(11) = 1; % rays intersect concave lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = mediumRefractiveIndex;
        opticalSystemOut = [opticalSystemOut; lensLine];
        
    case 'cameraToEye'
        
        % Add the front spectacle surface to the optical system.
        SlensFront = quadric.scale(quadric.unitSphere,[frontCurvature frontCurvature frontCurvature]);
        SlensFront = quadric.translate(SlensFront,[frontCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensFront);
        lensLine(11) = -1; % rays intersect convex lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = p.Results.lensRefractiveIndex;
        opticalSystemOut = [opticalSystemOut; lensLine];
        
        % Add the back spectacle surface to the optical system.
        SlensBack = quadric.scale(quadric.unitSphere,[backCurvature backCurvature backCurvature]);
        SlensBack = quadric.translate(SlensBack,[backCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensBack);
        lensLine(11) = -1; % rays intersect convex lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = mediumRefractiveIndex;
        opticalSystemOut = [opticalSystemOut; lensLine];
        
    otherwise
        error('This is not a valid setting for systemDirection');
end

end % function - addSpectacleLens