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
%  'lensRefractiveIndex'  - Scalar. Refractive index of the lens material.
%                           The routine returnRefractiveIndex() provides
%                           the indices for several spectacle materials
%                           under visible (vis) and near infra-red (nir)
%                           imaging domains.
%  'minimumLensThickness' - Scalar. The minimum lens thickness in mm.
%  'cornealRotation'      - 1x3 vector. If the cornea has been rotated out
%                           of alignment with the optical axis, the calling
%                           function needs to pass this value here. It
%                           would be better if I derived the rotation
%                           angles from the quadric itself, but this turns
%                           out to be complicated.
%
%
% Outputs:
%   opticalSystemOut      - An (m+1)x19 matrix, corresponding to the
%                           opticalSystemIn with the addition of the
%                           contact lens.
%
% Examples:
%{
    % Confirm that lenses of appropriate power are being created
    lensDiopters = 3;
    eye = modelEyeParameters('sphericalAmetropia',-2);
    opticalSystemIn = assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[]);
    eyePower = calcOpticalPower(opticalSystemIn);
    opticalSystemOut = addContactLens(opticalSystemIn,lensDiopters,'cornealRotation',eye.cornea.rotation);
    eyePowerWithLens = calcOpticalPower(opticalSystemOut);
    assert(abs(eyePowerWithLens - (eyePower + lensDiopters))<0.01);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('lensRefractionDiopters',@isnumeric);

% Optional
p.addParameter('lensRefractiveIndex',returnRefractiveIndex( 'hydrogel', 'NIR' ),@isnumeric);
p.addParameter('minimumLensThickness',0.05,@isnumeric);
p.addParameter('cornealRotation',[0, 0, 0],@isnumeric);

% parse
p.parse(opticalSystemIn, lensRefractionDiopters, varargin{:})


% Obtain the systemDirection for the passed optical system
if ~isempty(opticalSystemIn) && ~contains(calcSystemDirection(opticalSystemIn),'eyeToCamera')
    error('addContactLens:invalidSystemDirection','Lenses are only added to an optical system in the eyeToCamera direction')
end

% Detect the special case of lensRefractionDiopters == 0 and return the
% unmodified optical system
if lensRefractionDiopters==0
    opticalSystemOut = opticalSystemIn;
    return
end


%% Setup fixed lens paramters;

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.=
mediumRefractiveIndex = opticalSystemIn(end,end);

% Obtain the properties of the tear film from the optical system. This is
% the last surface. The index of refraction of the tear film is stored with
% the next-to-last surface
tearFilm = opticalSystemIn(end,:);
tearRefractiveIndex = opticalSystemIn(end-1,19);
tearS = tearFilm(1:10);
backCenter = quadric.center(tearS);
backCenter = backCenter(1);
backRadii = quadric.radii(tearS);
tearThickness = backRadii(1)+backCenter(1);

% The desired optical system will have its refractive power plus the called
% for lens refraction.
targetDiopters = calcOpticalPower(opticalSystemIn) + lensRefractionDiopters;


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
frontCenter = @(x) tearThickness + x(1) + x(2);

% Return the optical system for the candidate lens
mySystem = @(x) ...
    assembleLensSystem(opticalSystemIn, ...
    p.Results.lensRefractiveIndex, mediumRefractiveIndex, tearRefractiveIndex, ...
    x(1), frontCenter(x), ...
    -2, tearThickness, p.Results.cornealRotation);

% Calculate the power of the lens defined by the x parameters.
myDiopters = @(x) calcOpticalPower(mySystem(x));

% Define an objective which is the difference between the desired and
% measured optical power of the lens
myObj = @(x) ...
    abs(targetDiopters - myDiopters(x));

% Specify a non-linear constraint that requires that the lens has a
% positive radial thickness for a field of view that extends 63 degrees
% on either side of fixation. This is relevant only for plus lenses
myConstraint = @(x) checkLensShape(mySystem(x));

% Obtain an initial guess for the radius of curvature of the front
% surface using the thin lens approximation
frontCurvatureX0 = -backRadii(1);
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
    % upper bound on thickness. We allow the lower bound on thickness to go
    % below the "minimum", as otherwise we can't grind lenses that will
    % handle small (<1) positive corrections.
    x0 = [frontCurvatureX0*0.75 thicknessX0*2];
    lb = [frontCurvatureX0,p.Results.minimumLensThickness/4];
    ub = [frontCurvatureX0/4,thicknessX0*10];
    [x, fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,myConstraint,options);
else
    % This is a "minus" lens. Remove the non-linear shape constraint. Pin
    % the thickness to the minimum specified value.
    x0 = [frontCurvatureX0 thicknessX0];
    lb = [frontCurvatureX0*0.75,p.Results.minimumLensThickness];
    ub = [frontCurvatureX0/4,p.Results.minimumLensThickness];
    [x, fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,[],options);
end

if fVal > 0.01
    warning('addSpectacleLens:badGrind','Lens does not match requested optical power within tolerance. Perhaps hit a local minimum.');
end

% Restore the warning state
warning(warningState);


%% Assemble lens values
% In some cases overwrite the anonymous functions with the solution values
frontCurvature = x(1);
frontCenter = frontCenter(x);
[~,~,intersectHeight] = checkLensShape(mySystem(x));


%% Add the lens
opticalSystemOut =  assembleLensSystem(opticalSystemIn, p.Results.lensRefractiveIndex, mediumRefractiveIndex, tearRefractiveIndex, frontCurvature, frontCenter, intersectHeight, tearThickness, p.Results.cornealRotation);


end % function - addContactLens

%% LOCAL FUNCTIONS


function [c,ceq, intersectHeight] = checkLensShape(opticalSystem)

% Obtain the systemDirection
systemDirection = calcSystemDirection(opticalSystem);

% Use the systemDirection information to identify the "front" and "back"
% lenses from the perspective of the eye. This routine requires that the
% optical system is oriented "eyeToCamera"
switch systemDirection
    case 'eyeToCamera'
        % This is the front surface of the lens. The last surface is the
        % tear film
        Sfront = opticalSystem(end-1,1:10);
        % This is the front surface of the cornea.
        Sback = opticalSystem(end-3,1:10);
    case 'cameraToEye'
        error('addSpectacleLens:systemDirection','The routine only operates on optical systems oriented eyeToCamera');
    otherwise
        % We have an invalid optical system. This can occur if the
        % curvature of the lens becomes extreme. Return high constraint
        % vals to drive fmincon away from this part of the parameter space.
        c = -realmax;
        ceq = -realmax;
        return
end

% Calculate the distance from the corneal vertex to the back and front
% surface of the lens along a 63 degree viewing angle from the iris center.
R = quadric.normalizeRay(quadric.anglesToRay([-3.9;0;0], 63, 0 ));
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

% Return the position along the optical axis for the edge of the lens.
% This will be used subsequently to assemble a bounding box.
intersectHeight = Xfront(1);

end


function opticalSystemOut = assembleLensSystem(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, tearRefractiveIndex, frontCurvature, frontCenter, intersectHeight, tearThickness, cornealRotation)
% Assembles and returns an optical system matrix given input

% We are always operating in the 'eyeToCamera' system direction

% Obtain the radii of the last surface of the opticalSystemIn
backRadii = quadric.radii(opticalSystemIn(end,1:10));

% Create radii for the contact lens that account for the astigmatic
% ellipsoid form of the cornea
meanCurvDelta = mean(backRadii(2:3)+frontCurvature);
frontRadii = [-backRadii(1), -backRadii(2)-meanCurvDelta,  -backRadii(3)-meanCurvDelta];

% The opticalSystemIn ends with the outer surface of the tear film, from
% which the ray emerges into the refractive index of the medium. When a
% contact lens is added, the ray now emerges from the tear pool into the
% refractive index of the contact lens. We make that change here.
opticalSystemOut = opticalSystemIn;
opticalSystemOut(end,end) = lensRefractiveIndex;

% Define a bounding box for the front surface of the lens
boundingBoxLens = [intersectHeight frontCenter-frontCurvature -10 10 -10 10];

% Add the front contact lens surface to the optical system
SlensFront = quadric.scale(quadric.unitSphere,frontRadii);
SlensFront = quadric.rotate(SlensFront,cornealRotation);
SlensFront = quadric.translate(SlensFront,[frontRadii(1)+(frontCenter-frontCurvature) 0 0]);
lensLine = nan(1,19);
lensLine(1:10) = quadric.matrixToVec(SlensFront);
lensLine(11) = 1; % rays intersect concave lens surface
lensLine(12:17) = boundingBoxLens;
lensLine(18) = 1; % must intersect

% The ray emerges from the front contact lens surface into the
% refractive index of the tear film
lensLine(end) = tearRefractiveIndex;
opticalSystemOut = [opticalSystemOut; lensLine];

% Define a bounding box for the front tear film
boundingBoxTears = [intersectHeight+tearThickness frontCenter-frontCurvature+tearThickness -10 10 -10 10];

% Add a tear film to the optical system
SlensFront = quadric.scale(quadric.unitSphere,frontRadii);
SlensFront = quadric.rotate(SlensFront,cornealRotation);
SlensFront = quadric.translate(SlensFront,[frontRadii(1)+(frontCenter-frontCurvature)+tearThickness 0 0]);
tearLine = nan(1,19);
tearLine(1:10) = quadric.matrixToVec(SlensFront);
tearLine(11) = 1; % rays intersect concave lens surface
tearLine(12:17) = boundingBoxTears;
tearLine(18) = 1; % must intersect

% The ray emerges from the tear film into the refractice index of
% the medium
tearLine(end) = mediumRefractiveIndex;
opticalSystemOut = [opticalSystemOut; tearLine];

end


