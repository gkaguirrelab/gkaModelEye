function [opticalSystemOut, p] = addContactLens(opticalSystemIn, lensRefractionDiopters, varargin)
% Add a contact lens to a passed optical system
%
% Syntax:
%  [opticalSystemOut, p] = addContactLens(opticalSystemIn, lensRefractionDiopters)
%
% Description:
%	This routine adds a contact lens to an optical system with the
%	refractive power specified in the passed variable.
%
% Inputs:
%   opticalSystem         - An mx19 matrix, where m is set by the key value
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
% Optional key/value pairs:
%  'lensRefractiveIndex'  - Scalar. Refractive index of the lens material.
%                           The routine returnRefractiveIndex() provides
%                           the indices for several spectacle materials
%                           under visible (vis) and near infra-red (nir)
%                           imaging domains.
%  'minimumLensThickness' - Scalar. The minimum lens in mm.
%
%
% Outputs:
%   opticalSystemOut      - An (m+1)x19 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           contact lens.
%
% Examples:
%{
    lensDiopters = -3;
    eye = modelEyeParameters();
    opticalSystem=assembleOpticalSystem(eye,'surfaceSetName','stopToCamera','opticalSystemNumRows',[]);
    opticalSystem = addContactLens(opticalSystem,lensDiopters);
%}
%{
    %% Replicate calculation of contact lens curvature
    %% CURRENTLY NOT IMPLEMENTED
    % WJ Benajamin provudes an example calculation for a contact lens
    % that provides -10D power. We confirm here that our routine provides
    % the same solution.
    %   Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact 
    %   lens practice. Lippincott Williams & Wilkins, 2005. Chapter 7A, 
    %   "Optical phenomena of contact lenses", WJ Benjamin. p130
    % opticalSystemIn = [nan nan nan nan 1.3760; -7.8  -7.8  -7.8  -7.8  1.0];
    % The curvature of the front surface of the contact lens should be
    % -9.56. We obtain a slightly lower value (9.5468) as we do not model
    % the effect of the pre-lens tear film.
    % opticalSystemOut = addContactLens(opticalSystemIn, -10, 'lensRefractiveIndex', 1.43 );
    % assert(abs(opticalSystemOut(end,2) - -9.56)<0.1);
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

% Detect the special case of lensRefractionDiopters == 0 and return the
% unmodified optical system
if lensRefractionDiopters==0
    opticalSystemOut = opticalSystemIn;
    return
end

% Obtain the systemDirection for the passed optical system
systemDirection = calcSystemDirection(opticalSystemIn);
if ~strcmp(systemDirection,'eyeToCamera')
    error('Not set up yet')
end

%% Setup fixed lens paramters
% Distribute the parameters into variables
lensRefractiveIndex = p.Results.lensRefractiveIndex;

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.=
mediumRefractiveIndex = opticalSystemIn(end,end);

% Obtain the properties of the tear film from the optical system. This is
% the last surface. The index of refraction of the tear film is stored with
% the next-to-last surface
tearFilm = opticalSystemIn(end,:);
tearIndex = opticalSystemIn(end-1,19);
tearS = tearFilm(1:10);
backCenter = quadric.center(tearS);
backCenter = backCenter(1);
backCenter = quadric.center(tearS);
backCenter = backCenter(1);
backRadii = quadric.radii(tearS);
tearThickness = backRadii(1)+backCenter(1);

%% Search for parameters of an ophthalmic lens
% This is "convex-concave" lens with two surfaces. The back surface is
% equivalent to the tear film. Therefore, we need only find the front
% curvature and the position of the front surface. This search attempts to
% create a lens of the specified optical power, while satisfying a
% non-linear constraint upon the shape of the lens that varies depending
% upon whether the lens has positive or negative power.
%
% A plus lens corrects hyperopia. It is thickest along the optical axis. We
% constrain the shape of the lens so that it has no-zero thickness out at
% its far edge.
%
% A minus lens corrects myopia. It has a flatter front surface and a more
% curved back surface. It will be thinnest at the center of the lens on the
% optical axis. We constrain the shape of the lens so that the minimum
% thickness of the lens does not fall below a specified value.
%


%% Define anonymous functions
% Centers of the lenses, dependent upon lens thickness and back
% curvature
frontCenter = @(thickness) tearThickness + frontCurvature + thickness;

% Return the optical system for the candidate lens
myLens = @(x) ...
    assembleLensSystem(lensSystem, p.Results.systemDirection, ...
    p.Results.lensRefractiveIndex, mediumRefractiveIndex, ...
    x(1), backCenter(x(1)), ...
    frontCurvature, frontCenter(x(2)), ...
    [lensVertexDistance/2, lensVertexDistance/2], lensVertexDistance);

% Calculate the power of the lens defined by the x parameters
myDiopters = @(x) calcDiopters(myLens(x));

% Define an objective which is the difference between the desired and
% measured optical power of the lens
myObj = @(x) ...
    abs(lensRefractionDiopters - myDiopters(x));

% Specify a non-linear constraint that requires that the lens has a
% positive radial thickness for a field of view that extends 63 degrees
% on either side of fixation. This is relevant only for plus lenses
myConstraint = @(x) checkLensShape(myLens(x));

% Obtain an initial guess for the radius of curvature of the back
% surface using the thin lens approximation
backCurvatureX0 = sign(lensRefractionDiopters) * ...
    ((mediumRefractiveIndex-lensRefractiveIndex)/(lensRefractionDiopters - frontDiopters))*1000;
thicknessX0 = p.Results.minimumLensThickness;

% define some search options
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% Handle warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');


%% Perform the search
if sign(lensRefractionDiopters)==1
    % This is a "plus" lens. Apply the shape constraint but do not place an
    % upper bound on thickness.
    x0 = [backCurvatureX0 thicknessX0*2];
    lb = [-inf,p.Results.minimumLensThickness];
    ub = [inf,inf];
    [x, fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,myConstraint,options);
else
    % This is a "minus" lens. Remove the non-linear shape constraint. Pin
    % the thickness to the minimum specified value.
    x0 = [backCurvatureX0 thicknessX0];
    lb = [-inf,p.Results.minimumLensThickness];
    ub = [inf,p.Results.minimumLensThickness];
    [x, fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,[],options);    
end

if fVal > 0.01
    warning('addSpectacleLens:badGrind','Lens does not match requested optical power within tolerance');
end

% Restore the warning state
warning(warningState);


%% Assemble lens values
% In some cases overwrite the anonymous functions with the solution values
backCurvature = x(1);
backCenter = backCenter(backCurvature);
thickness = x(2);
frontCenter = frontCenter(thickness);
[~,~,intersectHeights] = checkLensShape(myLens(x));


%% Add the lens
opticalSystemOut = assembleLensSystem(opticalSystemIn, p.Results.systemDirection, p.Results.lensRefractiveIndex, mediumRefractiveIndex, backCurvature, backCenter, frontCurvature, frontCenter, intersectHeights, lensVertexDistance);

end % function - addSpectacleLens


%% LOCAL FUNCTIONS


function [c,ceq, intersectHeights] = checkLensShape(opticalSystem)

% Obtain the systemDirection
systemDirection = calcSystemDirection(opticalSystem);

% Use the systemDirection information to identify the "front" and "back"
% lenses from the perspective of the eye
switch systemDirection
    case 'eyeToCamera'
        Sfront = opticalSystem(3,1:10);
        Sback = opticalSystem(2,1:10);
    case 'cameraToEye'
        Sfront = opticalSystem(2,1:10);
        Sback = opticalSystem(3,1:10);
end

% Calculate the distance from the corneal vertex to the back and front
% surface of the lens along a 63 degree viewing angle.
R = quadric.normalizeRay(quadric.anglesToRay([0;0;0], 63, 0 ));
side = 1; % Our lenses are all concave w.r.t. a ray arising from eye
Xback = quadric.intersectRay(Sback,R,side);
Xfront = quadric.intersectRay(Sfront,R,side);
Dback = sqrt(sum(Xback.^2));
Dfront = sqrt(sum(Xfront.^2));

% This is the thickness of the lens at its edge. It is possible for this
% value to be negative for some values of lens curvature.
thicknessAtEdge = Dfront - Dback;

% The non-linear constraint values. fmincon will search for a solution in
% which c<=0 and c==0. That is, we would like the lens thickness at the
% edge to be close to zero and not negative.
c = -thicknessAtEdge;
ceq = c;

% Return the position along the optical axis for the edge of each surface.
% This will be used subsequently to assemble a bounding box.
intersectHeights = [Xback(1) Xfront(1)];

end


function opticalSystemOut = assembleLensSystem(opticalSystemIn, systemDirection, lensRefractiveIndex, mediumRefractiveIndex, backCurvature, backCenter, frontCurvature, frontCenter, intersectHeights, lensVertexDistance)
% Assembles and returns an optical system matrix given input

switch systemDirection
    case 'eyeToCamera'
        % Define a bounding box for the back surface
        boundingBoxLens = [intersectHeights(1) frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

        % Add the back spectacle surface to the optical system.        
        SlensBack = quadric.scale(quadric.unitSphere,[backCurvature backCurvature backCurvature]);
        SlensBack = quadric.translate(SlensBack,[backCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensBack);
        lensLine(11) = 1; % rays intersect concave lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = lensRefractiveIndex;
        opticalSystemOut = [opticalSystemIn; lensLine];

        % Define a bounding box for the front surface
        boundingBoxLens = [intersectHeights(2) frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

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
        % Define a bounding box for the front surface
        boundingBoxLens = [intersectHeights(2) frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

        % Add the front spectacle surface to the optical system.
        SlensFront = quadric.scale(quadric.unitSphere,[frontCurvature frontCurvature frontCurvature]);
        SlensFront = quadric.translate(SlensFront,[frontCenter 0 0]);
        lensLine = nan(1,19);
        lensLine(1:10) = quadric.matrixToVec(SlensFront);
        lensLine(11) = -1; % rays intersect convex lens surface
        lensLine(12:17) = boundingBoxLens;
        lensLine(18) = 1; % must intersect
        lensLine(end) = lensRefractiveIndex;
        opticalSystemOut = [opticalSystemIn; lensLine];
        
        % Define a bounding box for the back surface
        boundingBoxLens = [intersectHeights(1) frontCenter-frontCurvature -lensVertexDistance*2 lensVertexDistance*2 -lensVertexDistance*2 lensVertexDistance*2];

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


    
end % function - addContactLens