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
%  'backSurfaceRadii'     - 1x3 vector or empty. The radii of the tear film
%                           of the cornea, and thus the back surface of the
%                           contact lens. If left empty, these values are
%                           derived from the last row of the optical
%                           system.
%  'tearThickness'        - Scalar. The thickness of the tear film. If left
%                           empty, this value is derived from the last rows
%                           of the optical system
%  'contactLensDepth'     - Scalar. The axial cup depth (in mm) of the
%                           contact lens. Defaults to 2mm.
%  'cornealRotation'      - 1x3 vector. If the cornea has been rotated out
%                           of alignment with the optical axis, the calling
%                           function needs to pass this value here. It
%                           would be better if I derived the rotation
%                           angles from the quadric itself, but this turns
%                           out to be complicated.
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
    eye = modelEyeParameters('sphericalAmetropia',3);
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
p.addParameter('backSurfaceRadii',[],@isnumeric);
p.addParameter('tearThickness',[],@isnumeric);
p.addParameter('contactLensDepth',2,@isnumeric);
p.addParameter('contactLensViewAngle',63,@isnumeric);
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


%% Setup fixed lens paramters

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.
mediumRefractiveIndex = opticalSystemIn(end,end);

% Obtain the properties of the tear film from the optical system. This is
% the last surface. The index of refraction of the tear film is stored with
% the next-to-last surface
tearFilm = opticalSystemIn(end,:);
tearRefractiveIndex = opticalSystemIn(end-1,19);
tearS = tearFilm(1:10);
backCenter = quadric.center(tearS);
backCenter = backCenter(1);

% Derive the radii of the backSurface of the contact lens from the tear
% film, or the passed value.
if isempty(p.Results.backSurfaceRadii)
    backSurfaceRadii = quadric.radii(tearS);
else
    backSurfaceRadii = p.Results.frontSurfaceRadii;
end

if isempty(p.Results.tearThickness)
    tearThickness = backSurfaceRadii(1)+backCenter(1);
else
    tearThickness = p.Results.tearThickness;
end

% The desired optical system will have its refractive power plus the
% requested lens refraction.
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
% constrain the shape of the lens so that it has non-zero thickness at the
% far edge.
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

% Return the optical system for the candidate lens. We pass the negative of
% the contactLensDepth, as this defines the posterior aspect of the
% bounding box for the contact lens.
mySystem = @(x) ...
    assembleLensSystem(opticalSystemIn, ...
    p.Results.lensRefractiveIndex, mediumRefractiveIndex, ...
    tearRefractiveIndex, x(1), frontCenter(x), -p.Results.contactLensDepth, ...
    tearThickness, p.Results.cornealRotation, backSurfaceRadii);

% Define an objective which is the difference between the desired and
% measured optical power of the lens
myObj = @(x) ...
    abs(targetDiopters - myDiopters(x,mySystem));

% Specify a non-linear constraint that requires that the lens has a
% positive radial thickness for a field of view that extends 63 degrees
% on either side of fixation. This is relevant only for plus lenses
myConstraint = @(x) checkLensShape(mySystem(x),p.Results.contactLensViewAngle);

% Obtain an initial guess for the radius of curvature of the front
% surface using the thin lens approximation
frontCurvatureX0 = -backSurfaceRadii(1);
thicknessX0 = p.Results.minimumLensThickness;

% Define some search options
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
    x0 = [frontCurvatureX0*0.75, max([p.Results.minimumLensThickness, 0.042*lensRefractionDiopters])];
    lb = [frontCurvatureX0, p.Results.minimumLensThickness/4];
    ub = [frontCurvatureX0/4, thicknessX0*10];
    [x, fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,myConstraint,options);
else
    % This is a "minus" lens. Remove the non-linear shape constraint. Pin
    % the thickness to the minimum specified value.
    x0 = [frontCurvatureX0*0.75, thicknessX0];
    lb = [frontCurvatureX0, p.Results.minimumLensThickness];
    ub = [frontCurvatureX0/4, p.Results.minimumLensThickness];
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
[~,~,intersectHeight] = checkLensShape(mySystem(x),p.Results.contactLensViewAngle);


%% Add the lens
opticalSystemOut =  assembleLensSystem(opticalSystemIn, p.Results.lensRefractiveIndex, mediumRefractiveIndex, tearRefractiveIndex, frontCurvature, frontCenter, intersectHeight, tearThickness, p.Results.cornealRotation, backSurfaceRadii);


end % function - addContactLens


%% LOCAL FUNCTIONS

function diopters = myDiopters(x, mySystem)
% Calculate the power of the lens defined by the x parameters.

% mySystem is an anonymous function that provides the opticalSystem given
% the x parameters of the contact lens.
opticalSystem = mySystem(x);

% I had been encountering occasional "invalid system direction" errors in
% producing some contact lenses. Placing this calculation in a try-catch
% block to get some diagnostic information for this event. I think I solved
% the cause (it was due to corneal rotation not being properly modeled) but
% I have left this code here just in case the problem should return.
try
    diopters = calcOpticalPower(opticalSystem);
catch
    % Make the diopters something arbitrarily large, so that fmincon will
    % avoid whatever x parameters produced this error
    diopters = 1e6;
    % Print some diagnostic text to help figure out why we are getting
    % these bad traces
    fprintf('\n************************************\n')
    fprintf('addContactLens produced an invalid optical system\n')
    fprintf('Lens front curvature: %2.4f \n',x(1));
    fprintf('Lens thickness: %2.4f \n',x(2));
    fprintf('Optical system matrix:\n');
    opticalSystem
    fprintf('************************************\n\n')
end

end


function [c,ceq, intersectHeight] = checkLensShape(opticalSystem,contactLensViewAngle)

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
        % This is the tear film on the front surface of the cornea.
        Sback = opticalSystem(end-2,1:10);
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
% surface of the lens along the passed viewing angle from the iris center.
% We make this calculation in four different directions
thicknessAtEdge = [];
Xfront = [];
horiz = [-contactLensViewAngle contactLensViewAngle]; vert = [-contactLensViewAngle contactLensViewAngle];
for hh = 1:length(horiz)
    for vv = 1:length(vert)
        R = quadric.normalizeRay(quadric.anglesToRay([-3.9;0;0], horiz(hh), vert(hh) ));
        side = 1; % Our lenses are all concave w.r.t. a ray arising from eye
        Xback = quadric.intersectRayQuadric(Sback,R,side);
        Xfront(end+1,:) = quadric.intersectRayQuadric(Sfront,R,side);
        Dback = sqrt(sum(Xback.^2));
        Dfront = sqrt(sum(Xfront(end,:).^2));
        
        % This is the thickness of the lens at its edge. It is possible for
        % this value to be negative for some values of lens curvature.
        thicknessAtEdge(end+1) = Dfront - Dback;
    end
end

% The non-linear constraint values. fmincon will search for a solution in
% which c<=0 and ceq==0.

% The constraint is violated for any thickness that is negative
if any(thicknessAtEdge<0)
    c = -min(thicknessAtEdge);
else
    c = 0;
end

% The constraint is violated if the thinnest edge of the lens is greater
% than zero thickness
if any(thicknessAtEdge>0)
    ceq = min(thicknessAtEdge(thicknessAtEdge>0));
else
    ceq = 0;
end

% Return the position along the optical axis for the edge of the lens.
% This will be used subsequently to assemble a bounding box.
intersectHeight = min(Xfront(:,1));

end


function opticalSystemOut = assembleLensSystem(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, tearRefractiveIndex, frontCurvature, frontCenter, intersectHeight, tearThickness, cornealRotation, backRadii)
% Assembles and returns an optical system matrix given input. We are always
% operating in the 'eyeToCamera' system direction

% Create radii for the contact lens that account for the astigmatic
% ellipsoid form of the cornea
meanCurvDelta = mean(backRadii(2:3)+frontCurvature);
frontRadii = [-backRadii(1), -backRadii(2)-meanCurvDelta, -backRadii(3)-meanCurvDelta];

% The opticalSystemIn ends with the outer surface of the tear film, from
% which the ray emerges into the refractive index of the medium. When a
% contact lens is added, the ray now emerges from the tear pool into the
% refractive index of the contact lens. We make that change here.
opticalSystemOut = opticalSystemIn;
opticalSystemOut(end,end) = lensRefractiveIndex;

% Create the front contact lens surface
SlensFront = quadric.scale(quadric.unitSphere,frontRadii);
SlensFront = quadric.rotate(SlensFront,cornealRotation);
SlensFront = quadric.translate(SlensFront,[frontRadii(1)+(frontCenter-frontCurvature) 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( SlensFront );

% Define a bounding box for the front surface of the lens
boundingBoxLens = [intersectHeight X(1) -10 10 -10 10];

% Add this surface to the optical system
lensLine = nan(1,19);
lensLine(1:10) = quadric.matrixToVec(SlensFront);
lensLine(11) = 1; % rays intersect concave lens surface
lensLine(12:17) = boundingBoxLens;
lensLine(18) = 1; % must intersect

% The ray emerges from the front contact lens surface into the refractive
% index of the tear film
lensLine(end) = tearRefractiveIndex;
opticalSystemOut = [opticalSystemOut; lensLine];

% Add a tear film to the optical system
SlensTear = quadric.scale(quadric.unitSphere,frontRadii);
SlensTear = quadric.rotate(SlensTear,cornealRotation);
SlensTear = quadric.translate(SlensTear,[frontRadii(1)+(frontCenter-frontCurvature)+tearThickness 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( SlensTear );

% Define a bounding box for the front tear film
boundingBoxTears = [intersectHeight+tearThickness X(1) -10 10 -10 10];

% Add this surface to the optical system
tearLine = nan(1,19);
tearLine(1:10) = quadric.matrixToVec(SlensTear);
tearLine(11) = 1; % rays intersect concave lens surface
tearLine(12:17) = boundingBoxTears;
tearLine(18) = 1; % must intersect

% The ray emerges from the tear film into the refractive index of the
% medium
tearLine(end) = mediumRefractiveIndex;
opticalSystemOut = [opticalSystemOut; tearLine];

end


