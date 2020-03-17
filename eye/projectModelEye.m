function [pupilEllipseOnImagePlane, glintCoord, imagePoints, worldPoints, headPoints, eyePoints, pointLabels, targetIntersectError, pupilFitError] = projectModelEye(eyePose, sceneGeometry, varargin)
% Obtain the parameters of the entrance pupil ellipse on the image plane
%
% Syntax:
%  [pupilEllipseOnImagePlane, glintCoord, imagePoints, worldPoints, headPoints, eyePoints, pointLabels, targetIntersectError, pupilFitError] = projectModelEye(eyePose, sceneGeometry)
%
% Description:
%   Given the sceneGeometry, this routine simulates the aperture stop in a
%   rotated eye and returns the parameters of the ellipse (in transparent
%   format) fit to the entrance pupil in the image plane. The location of
%   any glint(s) are also provided.
%
%   The forward model is a perspective projection of an anatomically
%   accurate eye, with points positioned behind the cornea subject to
%   refractive displacement. The projection incorporates the intrinsic
%   properties of the camera, including any radial lens distortion.
%
%   The 'fullEyeModelFlag' key will render the entire eye, not just the
%   pupil.
%
% Notes:
%   Rotations - Eye rotation is given as azimuth, elevation, and torsion in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds to
%   the position of the eye when the optical axis of the eye is aligned
%   with the optical axis of the camera. Positive rotations correspond to
%   rightward / upward (+x, -y) translation of the pupil center in the
%   image. Torsion of zero corresponds to the torsion of the eye when it is
%   in primary position.
%
%   Units - Eye rotations are in units of degrees. However, the units of
%   theta in the transparent ellipse parameters are radians. This is in
%   part to help us keep the two units separate conceptually.
%
% Inputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'addPseudoTorsion'     - Logical. If set to true, the eyePose is
%                           adjusted to add "pseudo" torsion so that the
%                           eye movement obeys Listing's Law. More details
%                           here:   addPseudoTorsion.m
%  'fullEyeModelFlag'     - Logical. Determines if the full eye model will
%                           be created.
%  'nStopPerimPoints'     - Scalar. The number of points that are
%                           distributed around the stop ellipse. A minimum
%                           of 5 is required to uniquely specify the image
%                           ellipse, and 6 to obtain a meaningful
%                           pupilFitError.
%  'stopPerimPhase'       - Scalar. The phase (in radians) of the position
%                           of the stop perimeter points.
%  'replaceReflectedPoints' - Logical. Points that are subject to
%                           refraction may encounter total internal
%                           reflection and therefore not appear in the
%                           image. If set to true, if the point is on the
%                           stop or iris border, then the routine will
%                           search to find the closest point from the pupil
%                           or iris that would instead appear in the image.
%                           This supports the non-elliptical appearance of
%                           the iris (and to a lesser extent the pupil) at
%                           extreme viewing angles.
%  'borderSearchPrecision' - Scalar. When the search process defined by the
%                           replaceReflectedPoints flags is performed, this
%                           parameter defines the precision of the search.
%  'rayTraceErrorThreshold' - Scalar. Virtual image points that have
%                           a ray trace error above this threshold will be
%                           discarded.
%  'nIrisPerimPoints'     - Scalar. The number of points that are
%                           distributed around the iris circle.
%  'corneaMeshDensity'    - Scalar. The number of geodetic lines used to
%                           render the corneal ellipsoid. About 20 makes a
%                           nice image.
%  'retinaMeshDensity'    - Scalar. The number of geodetic lines used to
%                           render the retina ellipsoid. About 24 makes a
%                           nice image.
%  'pupilRayFunc','glintRayFunc' - Function handles. By default, these are 
%                           set to 'findPupilRayMex' and 'findGlintRayMex'.
%                           This option is provided so that the routine can
%                           be tested with the native MATLAB code.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
%   glintCoord            - A nx2 vector with the image coordinates of the
%                           n glints.
%   imagePoints           - An nx2 matrix that specifies the x, y location
%                           on the image plane for each of the eyeWorld
%                           points.
%   worldPoints           - An nx3 matrix of the coordinates of the
%                           points of the eye model in the world
%                           coordinate frame. If fullEyeModel is set to
%                           false, then there will only be points for the
%                           pupil perimeter. If fullEyeModel is true, then
%                           the entire model of ~1000 points will be
%                           returned.
%   headPoints            - An nx3 matrix of the coordinates of the
%                           points of the eye model in the headWorld
%                           coordinate frame. This is the same as the
%                           eyePoints coordinate frame after eye rotation.
%   eyePoints             - An nx3 matrix of the coordinates of the
%                           points of the eye model in the eyeWorld
%                           coordinate frame.
%   pointLabels           - An nx1 cell array that identifies each of the
%                           points, with values such as: {'stopCenter',
%                           'irisCenter', 'aziRotationCenter',
%                           'eleRotationCenter', 'retina',
%                           'irisPerimeter', 'pupilPerimeter',
%                           'cornea','cornealApex'}.
%   targetIntersectError -  A nx1 vector that contains the distance (in mm)
%                           between the pinhole aperture of the camera and
%                           the intersection of a ray on the camera plane.
%                           The value is nan for points not subject to
%                           refraction by the cornea. All values will be
%                           nan if sceneGeometry.refraction is empty.
%   pupilFitError         - The RMSE distance (in pixels) of the pupil
%                           perimeter points to the pupil ellipse. This
%                           value indicates the degree to which the shape
%                           of the entrance pupil departs from an ellipse.
%
% Examples:
%{
    %% Basic forward projection example
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with the azimuth, elevation, torsion, and stop radius
    eyePose = [-10 5 0 3];
    % Obtain the pupil ellipse parameters in transparent format for the
    % default sceneGeometry
    pupilEllipseOnImagePlane = projectModelEye(eyePose,sceneGeometry);
    % Test against cached result
    pupilEllipseOnImagePlaneCached = [ 0.027801112018072   0.022387973426570   1.554201447325573   0.000023117837874   0.000191404712119 ].*1e4;
    assert(max(abs(pupilEllipseOnImagePlane -  pupilEllipseOnImagePlaneCached)) < 1e-4)
%}
%{
    %% Test the accuracy of the ellipse fit to the pupil boundary
    sceneGeometry=createSceneGeometry();
    aziVals = -60:5:60;
    pupilFitError = [];
    for aa = 1:length(aziVals)
        eyePose = [aziVals(aa) -3 0 3];
        [pupilEllipseOnImagePlane, ~, ~, ~, ~, ~, ~, ~, pupilFitError(aa)] = projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16);
    end
    figure
    plot(aziVals,pupilFitError,'.r');
%}
%{
    %% Show the non-elliptical iris perimeter at extreme viewing angles
    sceneGeometry=createSceneGeometry();
    eyePose = [-65 0 0 3];
    [pupilEllipseOnImagePlane, ~, imagePoints, ~, ~, ~, pointLabels] = ...
        projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16, ...
        'replaceReflectedPoints',true, ...
        'nIrisPerimPoints',16,'fullEyeModelFlag', true);
    figure
    idx = strcmp(pointLabels,'pupilPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'.r');
    hold on
    addTransparentEllipseToFigure(pupilEllipseOnImagePlane);
    idx = strcmp(pointLabels,'irisPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'.b');
    [pupilEllipseOnImagePlane, ~, imagePoints, ~, ~, ~, pointLabels] = ...
        projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16, ...
        'replaceReflectedPoints',false, ...
        'nIrisPerimPoints',16,'fullEyeModelFlag', true);
    idx = strcmp(pointLabels,'pupilPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'xr');
    addTransparentEllipseToFigure(pupilEllipseOnImagePlane);
    idx = strcmp(pointLabels,'irisPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'xb');
    axis equal
%}
%{
    %% Calculate the time required for the forward projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Generate some randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*20, (rand(nPoses,1)-0.5)*10, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    fprintf('Time to compute forward projection model (average over %d projections):\n',nPoses);
    tic
    for pp = 1:nPoses
    	projectModelEye(eyePoses(pp,:),sceneGeometry);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing compiled ray tracing: %4.2f msecs.\n',msecPerModel);
    tic
    for pp = 1:nPoses
    	projectModelEye(eyePoses(pp,:),sceneGeometry,'pupilRayFunc',@findPupilRay);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing MATLAB ray tracing: %4.2f msecs.\n',msecPerModel);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eyePose',@(x)(isnumeric(x) && all(size(x)==[1 4])));
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('addPseudoTorsion',true,@islogical);
p.addParameter('fullEyeModelFlag',false,@islogical);
p.addParameter('nStopPerimPoints',6,@(x)(isscalar(x) && x>4));
p.addParameter('stopPerimPhase',0,@isscalar);
p.addParameter('replaceReflectedPoints',false,@islogical);
p.addParameter('borderSearchPrecision',0.01,@isscalar);
p.addParameter('rayTraceErrorThreshold',0.01,@isscalar);
p.addParameter('nIrisPerimPoints',5,@isscalar);
p.addParameter('corneaMeshDensity',23,@isscalar);
p.addParameter('retinaMeshDensity',30,@isscalar);
p.addParameter('pupilRayFunc',@findPupilRayMex,@(x)(isa(x,'function_handle') || isempty(x)));
p.addParameter('glintRayFunc',@findGlintRayMex,@(x)(isa(x,'function_handle') || isempty(x)));

% parse
p.parse(eyePose, sceneGeometry, varargin{:})


%% Extract fields from Results struct
% Accessing struct elements is relatively slow. Do this here so it is not
% done in a loop.
pupilRayFunc = p.Results.pupilRayFunc;
glintRayFunc = p.Results.glintRayFunc;
replaceReflectedPoints = p.Results.replaceReflectedPoints;
rayTraceErrorThreshold = p.Results.rayTraceErrorThreshold;
borderSearchPrecision = p.Results.borderSearchPrecision;


%% Apply pseudoTorsion if requested
if p.Results.addPseudoTorsion
    eyePose = addPseudoTorsion(eyePose,sceneGeometry.eye.rotationCenters.primaryPosition);
end


%% Prepare variables
stopRadius = eyePose(4);
glintCoord = [];

% Store the number of pupil perimeter points
nStopPerimPoints = p.Results.nStopPerimPoints;


%% Define an eye in eye coordinates
% This coordinate frame is in mm and has the dimensions (p1,p2,p3). The
% diagram is of a cartoon pupil, viewed directly from the front.
%
% Coordinate [0,0,0] corresponds to the apex (front surface) of the cornea,
% with the model eye having the property of the optical and pupil axes of
% the eye being aligned. The first dimension is depth, and has a negative
% value toward the back of the eye.
%
%                 |
%     ^         __|__
%  +  |        /     \
% p3  -  -----(   +   )-----
%  -  |        \_____/
%     v           |
%                 |
%
%           - <--p2--> +
%
% For the right eye, negative values on the p2 dimension are more temporal,
% and positive values are more nasal. Positive values of p3 are upward,
% and negative values are downward


%% Define points around the elliptical aperture stop

% The eccentricity of the pupil aperture is given by a stored function
stopEccenFunc = str2func(sceneGeometry.eye.stop.eccenFcnString);

% Determine the parameters of the ellipse that defines the aperture stop in
% the plane of the iris. The absolute value of stopEccenFunc gives the
% eccentricity. The theta of the stop switches from horizontal to vertical
% when the stop passes from a negative to positive eccentricity, passing
% through circular at an eccentricity of 0.
stopEccen = stopEccenFunc(stopRadius);
stopEllipse = [sceneGeometry.eye.stop.center(2) , ...
    sceneGeometry.eye.stop.center(3), ...
    pi*stopRadius^2, ...
    abs(stopEccen),...
    sceneGeometry.eye.stop.thetas(1+(stopEccen>0))];

% Obtain the points on the perimeter of the ellipse
[p2p, p3p] = ellipsePerimeterPoints( stopEllipse, nStopPerimPoints, p.Results.stopPerimPhase );

% Place these points into the eyeWorld coordinates. Optionally create
% separate front and back stop perimeters to model iris thickness.
if sceneGeometry.eye.iris.thickness~=0
    stopPoints = zeros(nStopPerimPoints*2,3);
    stopPoints(1:nStopPerimPoints*2,3) = [p3p; p3p];
    stopPoints(1:nStopPerimPoints*2,2) = [p2p; p2p];
    stopPoints(1:nStopPerimPoints,1) = sceneGeometry.eye.stop.center(1)+sceneGeometry.eye.iris.thickness/2;
    stopPoints(nStopPerimPoints+1:nStopPerimPoints*2,1) = sceneGeometry.eye.stop.center(1)-sceneGeometry.eye.iris.thickness/2;
    
    % Create labels for the stopPerimeter points
    tmpLabelsFront = repmat({'stopPerimeterFront'},nStopPerimPoints, 1);
    tmpLabelsBack = repmat({'stopPerimeterBack'},nStopPerimPoints, 1);
    pointLabels = [tmpLabelsFront; tmpLabelsBack];
else
    stopPoints = zeros(nStopPerimPoints,3);
    stopPoints(1:nStopPerimPoints,3) = p3p;
    stopPoints(1:nStopPerimPoints,2) = p2p;
    stopPoints(1:nStopPerimPoints,1) = sceneGeometry.eye.stop.center(1);
    
    % Create labels for the stopPerimeter points
    tmpLabels = repmat({'stopPerimeter'},nStopPerimPoints, 1);
    pointLabels = tmpLabels;
end
eyePoints = stopPoints;


%% Define full eye model
% If the fullEyeModel flag is set, then we will create a model of the
% posterior and anterior segments of the eye.
if p.Results.fullEyeModelFlag
    
    % Add points for the stop center, iris center, rotation centers,
    % and retinal landmarks if present
    eyePoints = [eyePoints; sceneGeometry.eye.stop.center];
    pointLabels = [pointLabels; 'stopCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.iris.center];
    pointLabels = [pointLabels; 'irisActualCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.azi];
    pointLabels = [pointLabels; 'aziRotationCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.ele];
    pointLabels = [pointLabels; 'eleRotationCenter'];
    eyePoints = [eyePoints; 0 0 0];
    pointLabels = [pointLabels; 'vertex'];
    if isfield(sceneGeometry.eye,'landmarks')
        if isfield(sceneGeometry.eye.landmarks,'fovea')
            eyePoints = [eyePoints; sceneGeometry.eye.landmarks.fovea.coords];
            pointLabels = [pointLabels; 'fovea'];
        end
        if isfield(sceneGeometry.eye.landmarks,'opticDisc')
            eyePoints = [eyePoints; sceneGeometry.eye.landmarks.opticDisc.coords];
            pointLabels = [pointLabels; 'opticDisc'];
        end
    end
    
    % Define points around the perimeter of the iris
    nIrisPerimPoints = p.Results.nIrisPerimPoints;
    perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
    irisPoints = zeros(nIrisPerimPoints,3);
    irisPoints(1:nIrisPerimPoints,3) = ...
        sin(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(3);
    irisPoints(1:nIrisPerimPoints,2) = ...
        cos(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(2);
    irisPoints(1:nIrisPerimPoints,1) = ...
        0 + sceneGeometry.eye.iris.center(1);
    
    % Add the points and labels
    eyePoints = [eyePoints; irisPoints];
    tmpLabels = cell(size(irisPoints,1), 1);
    tmpLabels(:) = {'irisActualPerimeter'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Create the corneal surface points
    corneaPoints = quadric.surfaceGrid(...
        sceneGeometry.eye.cornea.front.S,...
        sceneGeometry.eye.cornea.front.boundingBox,...
        p.Results.corneaMeshDensity, ...
        'parametricPolar');
    
    % Identify the index of the corneal apex
    [~,apexIdx]=max(corneaPoints(:,1));
    
    % Save the corneal apex coordinates
    cornealApex = corneaPoints(apexIdx,:);
    
    % Add the points and labels
    eyePoints = [eyePoints; corneaPoints];
    tmpLabels = cell(size(corneaPoints,1), 1);
    tmpLabels(:) = {'cornea'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Add an entry for the corneal apex
    eyePoints = [eyePoints; cornealApex];
    pointLabels = [pointLabels; 'cornealApex'];
    
    % Create the retina surface vertices
    retinaPoints = quadric.surfaceGrid(...
        sceneGeometry.eye.retina.S,...
        sceneGeometry.eye.retina.boundingBox,...
        p.Results.retinaMeshDensity, ...
        'ellipsoidalPolar');
    
    % Retain those points that are posterior to the iris plane, and have a
    % distance from the optical axis in the p2xp3 plane of greater than the
    % iris radius
    retainIdx = logical(...
        (retinaPoints(:,1) < sceneGeometry.eye.iris.center(1)) .* ...
        sqrt(retinaPoints(:,2).^2+retinaPoints(:,3).^2) > sceneGeometry.eye.iris.radius );
    if all(~retainIdx)
        error('projectModelEye:irisCenterPosition','The iris center is behind the center of the retina');
    end
    retinaPoints = retinaPoints(retainIdx,:);
    
    % Add the points and labels
    eyePoints = [eyePoints; retinaPoints];
    tmpLabels = cell(size(retinaPoints,1), 1);
    tmpLabels(:) = {'retina'};
    pointLabels = [pointLabels; tmpLabels];
    
end


%% Refract the eyeWorld points
% This steps accounts for the effect of corneal and corrective lens
% refraction upon the appearance of points from the eye.

% Create a nan vector for the ray trace errors
targetIntersectError = nan(size(eyePoints,1),1);

% Identify the eyePoints subject to refraction by the cornea
refractPointsIdx = find(...
    strcmp(pointLabels,'stopPerimeter')+...
    strcmp(pointLabels,'stopPerimeterFront')+...
    strcmp(pointLabels,'stopPerimeterBack')+...
    strcmp(pointLabels,'irisActualPerimeter'));

% Check if we have a refraction field and it is not empty
refractFlag = false;
if isfield(sceneGeometry,'refraction') && ~isempty(pupilRayFunc)
    if ~isempty(sceneGeometry.refraction)
        refractFlag = true;
    end
end

% Perform refraction
if refractFlag
    % Assemble the static args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
        sceneGeometry.eye.rotationCenters, ...
        sceneGeometry.refraction.stopToMedium.opticalSystem, ...
        sceneGeometry.refraction.mediumToCamera.opticalSystem};
    
    % Pre-allocate variables to hold the results
    virtualPoints = nan(length(refractPointsIdx),3);
    virtualIntersectError = inf(length(refractPointsIdx),1);
    virtualPointLabels = cell(length(refractPointsIdx),1);
    virtualPointLabels(:) = {''};
    
    % Loop through the eyePoints that are to be refracted
    for ii=1:length(refractPointsIdx)
        
        % Get this eye point
        eyePoint=eyePoints(refractPointsIdx(ii),:);
        
        % Perform the computation using the passed function handle.
        [virtualImageRay, ~, intersectError] = ...
            pupilRayFunc(eyePoint, eyePose, args{:});
        virtualPoint = virtualImageRay(1,:);
        
        % Check if the point has encountered total internal reflection or
        % is a bad ray trace
        retainPoint = true;        
        
        if isnan(intersectError) || (intersectError > rayTraceErrorThreshold)
            % The eyePoint did not yield a valid image point. We will
            % not retain the point unless we find a replacement.
            retainPoint = false;
            
            % If this eyePoint is on the stop border, search across
            % smaller stop radii to find a replacement point on the
            % boundary of the pupil perimeter that does make it
            % through.
            if replaceReflectedPoints                
                % Find the appropriate center target for the search
                switch pointLabels{refractPointsIdx(ii)}
                    case 'stopPerimeter'
                        centerTarget = sceneGeometry.eye.stop.center;
                    case 'stopPerimeterFront'
                        centerTarget = sceneGeometry.eye.stop.center(1)+sceneGeometry.eye.iris.thickness/2;
                    case 'stopPerimeterBack'
                        centerTarget = sceneGeometry.eye.stop.center(1)-sceneGeometry.eye.iris.thickness/2;
                    case 'irisActualPerimeter'
                        centerTarget = sceneGeometry.eye.iris.center;
                end
                % Initialize the searchScalar
                searchScalar = 1.0 - borderSearchPrecision;
                % Perform the search
                stillSearching = true;
                while stillSearching
                    % Move the eyePoint along the vector connecting the eye
                    % point to the stop or iris center
                    shiftedEyePoint = eyePoint - (eyePoint - centerTarget).*(1 - searchScalar);
                    % Subject the shifted point to refraction
                    [virtualImageRay, ~, intersectError] = ...
                        pupilRayFunc(shiftedEyePoint, eyePose, args{:});
                    virtualPoint = virtualImageRay(1,:);
                    % Update the search scalar
                    searchScalar = searchScalar - borderSearchPrecision;
                    if searchScalar <= 0
                        stillSearching = false;
                    end
                    % Test if the newly refracted point meets criterion
                    if intersectError < rayTraceErrorThreshold
                        stillSearching = false;
                        retainPoint = true;
                    end
                end
            else
                retainPoint = false;
            end
        end
        
        % If the ray-trace process yielded a point in the image, store the
        % result
        if retainPoint
            
            % Add the refracted point to the set
            virtualPoints(ii,:) = virtualPoint;
            virtualIntersectError(ii) = intersectError;
            
            % Create a label for the virtual image point
            newPointLabel = pointLabels{refractPointsIdx(ii)};
            newPointLabel = strrep(newPointLabel,'stop','pupil');
            newPointLabel = strrep(newPointLabel,'Actual','');
            virtualPointLabels(ii) = {newPointLabel};
            
        end
    end
    
    % Add the virtual image points to the entire set of eye points. First
    % find those points that had a bad ray trace
    badTraceIdx = any(isnan(virtualPoints),2);
    
    % Remove the bad ray trace points from both the original eyePoint
    % set and the virtual image set
    if any(badTraceIdx)
        virtualPoints(badTraceIdx,:) = [];
        virtualPointLabels(badTraceIdx) = [];
        virtualIntersectError(badTraceIdx) = [];
        
        eyePoints(refractPointsIdx(badTraceIdx),:) = [];
        pointLabels(refractPointsIdx(badTraceIdx)) = [];
        targetIntersectError(refractPointsIdx(badTraceIdx)) = [];
    end
    eyePoints = [eyePoints; virtualPoints];
    pointLabels = [pointLabels; virtualPointLabels];
    targetIntersectError = [targetIntersectError; virtualIntersectError];
    
else
    % If there is no refraction, then the pupil is simply the stop. Copy
    % these points over to their new names
    for ii=1:length(refractPointsIdx)
        
        % Get this eyeWorld point
        eyePoint=eyePoints(refractPointsIdx(ii),:);
        
        % Add the refracted point to the set
        eyePoints = [eyePoints; eyePoint];
        targetIntersectError = [targetIntersectError; 0];
        
        % Create a label for the virtual image point
        newPointLabel = pointLabels{refractPointsIdx(ii)};
        newPointLabel = strrep(newPointLabel,'stop','pupil');
        newPointLabel = strrep(newPointLabel,'Actual','');
        pointLabels = [pointLabels; {newPointLabel}];
    end
    
end


%% Calc glint
% Find the location of the glint in the eye world coordinate frame. The
% glint is the reflection of a light source from the tear film of the eye.
% The location of the glint in the image is subject to refraction by
% artificial lenses.
if isfield(sceneGeometry,'refraction') && ~isempty(glintRayFunc)
    if isfield(sceneGeometry.refraction,'glint')
        
        % The position of the light source in the world coordinate frame
        % that is the source of the glint
        glintSourceWorld = sceneGeometry.cameraPosition.translation + ...
            sceneGeometry.cameraPosition.glintSourceRelative;
        
        % How many light sources do we have?
        nGlints = size(glintSourceWorld,2);
        
        % Assemble the args
        args = {sceneGeometry.cameraPosition.translation, ...
            sceneGeometry.eye.rotationCenters, ...
            sceneGeometry.refraction.cameraToMedium.opticalSystem, ...
            sceneGeometry.refraction.glint.opticalSystem, ...
            sceneGeometry.refraction.mediumToCamera.opticalSystem};
        
        % Loop through the glints
        for gg = 1:nGlints
            
            % Perform the computation using the passed function handle
            [virtualImageRay, ~, intersectError] = ...
                glintRayFunc(glintSourceWorld(:, gg), eyePose, args{:});
            
            % If we have a good trace, add the glint point and label
            if intersectError < rayTraceErrorThreshold
                eyePoints = [eyePoints; virtualImageRay(1,:)];
                glintLabel = sprintf('glint_%02d',gg);
                pointLabels = [pointLabels; {glintLabel}];
            end
            
        end
    end
end


%% Apply eye rotation

% Omit the eye rotation centers themselves from rotation.
rotatePointsIdx = find(~contains(pointLabels,{'Rotation'}));

% Copy the eyePoints to the headPoints
headPoints = eyePoints;

% Loop through the points to be rotated. We pass the rotation matrix to
% avoid having to re-calculate this for the rotation of each point.
R = [];
for pp = 1:length(rotatePointsIdx)
    headPoints(rotatePointsIdx(pp),:) = rotateEyeCoord(...
        eyePoints(rotatePointsIdx(pp),:), ...
        eyePose, ...
        sceneGeometry.eye.rotationCenters, ...
        'forward', ...
        R);
end

% If we are projecting a full eye model, label as hidden those posterior
% segment points that are posterior to the most posterior of the centers of
% rotation of the eye, and thus would not be visible to the camera.
if p.Results.fullEyeModelFlag
    seenIdx = strcmp(pointLabels,'retina') .* (headPoints(:,1) >= min([sceneGeometry.eye.rotationCenters.azi(1) sceneGeometry.eye.rotationCenters.ele(1)]));
    seenIdx = logical(seenIdx + ~strcmp(pointLabels,'retina'));
    pointLabels(~seenIdx) = strcat(pointLabels(~seenIdx),'_hidden');
end


%% Switch axes to world coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (taken from Leszek Swirski), being
% viewed from above:
%
%    ^
%    |
%    |    .-.
% -Z |   |   | <- Head
%    +   `^u^'
% +Z |
%    |
%    |      W <- Camera    (As seen from above)
%    V
%
%     <-----+----->
%        -X   +X
%
% +X = right
% +Y = up
% +Z = front (towards the camera)
%
% The origin [0,0,0] corresponds to the front surface of the eye and the
% pupil center when the optical axis of the eye and the camera axis are
% aligned.

% Re-arrange the headPoints to transform to the world coordinate frame
worldPoints = headPoints(:,[2 3 1]);


%% Project the world coordinate points to the image plane
% This coordinate frame is in units of pixels, and has the dimensions
% [x, y]:
%
% [0,0]    x
%      +------->
%      |
%   y  |
%      |
%      v
%
% With x being left/right and y being up/down
%

% Return if the sceneGeometry is lacking a cameraPosition or
% cameraIntrinsic field
if ~isfield(sceneGeometry,'cameraPosition') || ~isfield(sceneGeometry,'cameraIntrinsic')
    imagePoints = [];
    pupilEllipseOnImagePlane=nan(1,5);
    pupilFitError = nan;
    return
end

% Create the camera position rotation matrix. This is the rotation matrix
% for the position of the camera with respect to the world coordinates.
% Only camera torsion is supported.
cameraRotationMatrix = ...
    [cosd(sceneGeometry.cameraPosition.torsion)	-sind(sceneGeometry.cameraPosition.torsion)	0; ...
    sind(sceneGeometry.cameraPosition.torsion)     cosd(sceneGeometry.cameraPosition.torsion)     0; ...
    0                                   0                                       1];

% We need to post-multiply the camera rotation matrix by a matrix that
% reverses the Y and Z axes. This is to allow the camera coordinates to
% also be right handed (as are the world coordinates), and to accomodate
% the reversed Y-axis (negative is up) of the image plane.
cameraRotationMatrix = cameraRotationMatrix * [1 0 0; 0 -1 0; 0 0 -1];

% Create the extrinsic matrix. The sceneGeometry structure specifies the
% translation and rotation parameters of camera in world coordinates (i.e.,
% the sceneGeometry structure specifies the camera pose). The extrinsic
% matrix is equal to [R t], where R is transpose of the
% cameraRotationMatrix, and t = -RC, where C is the camerra translation in
% world coordinates
cameraExtrinsicMatrix = ...
    [transpose(cameraRotationMatrix) -transpose(cameraRotationMatrix)*sceneGeometry.cameraPosition.translation];

% The projection matrix is the cameraInstrinsic matrix times the
% cameraExtrinsic matrix.
projectionMatrix = sceneGeometry.cameraIntrinsic.matrix * cameraExtrinsicMatrix;

% What is our total number of points to project?
nEyePoints = size(eyePoints,1);

% Project the world points to camera points, and then to the image plane.
% The world points have a column of ones added to support the
% multiplication with a combined rotation and translation matrix
cameraPoints=(projectionMatrix*[worldPoints, ones(nEyePoints,1)]')';
imagePointsPreDistortion=zeros(nEyePoints,2);
imagePointsPreDistortion(:,1) = ...
    cameraPoints(:,1)./cameraPoints(:,3);
imagePointsPreDistortion(:,2) = ...
    cameraPoints(:,2)./cameraPoints(:,3);


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the sensor optical center (aka
% principal point), and the coordinates are in world units. To apply this
% distortion to our image coordinate points, we subtract the optical
% center, and then divide by fx and fy from the intrinsic matrix.

% Implement as bsxfun to avoid implicit expansion
imagePointsNormalized = bsxfun(@rdivide, ...
    bsxfun(@minus,imagePointsPreDistortion,[sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]) , ...
    [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]);

% Distortion is proportional to distance from the center of projection on
% the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(1).*radialPosition.^2 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(2).*radialPosition.^4;

imagePointsNormalizedDistorted = imagePointsNormalized.*0;
imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

% Place the distorted points back into the imagePoints vector
imagePoints = bsxfun(@plus, ...
    bsxfun(@times,imagePointsNormalizedDistorted , [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]) ,...
    [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]);


%% Store the glintCoord
glintIdx = strncmp(pointLabels,'glint',5);
if any(glintIdx)
    glintCoord = imagePoints(glintIdx,:);
else
    glintCoord = [];
end


%% Obtain the pupil ellipse
% Proceed with fitting if we have a non-zero stop radius
if eyePose(4) > 0
    if sceneGeometry.eye.iris.thickness==0
        % The simple case of a zero-thickness aperture stop. Identify the
        % perimeter points.
        pupilPerimIdx = logical(strcmp(pointLabels,'pupilPerimeter'));
        % Fit the ellipse
        [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
    else
        % The more complicated case of a non-zero thickness aperture stop.
        % Identify the perimeter points for the front and back virtual
        % image of the aperture stop.
        pupilPerimIdxFront = logical(strcmp(pointLabels,'pupilPerimeterFront'));
        pupilPerimIdxBack = logical(strcmp(pointLabels,'pupilPerimeterBack'));
        % Obtain the ellipse fit to the front and back
        frontEllipse = pupilEllipseFit(imagePoints(pupilPerimIdxFront,:));
        backEllipse = pupilEllipseFit(imagePoints(pupilPerimIdxBack,:));
        % If either ellipse fit yielded nans, then the final ellipse is nan
        if any(isnan(frontEllipse)) || any(isnan(backEllipse))
            % If we are unable to fit an ellipse to the front or back pupil
            % perimeter, then exit with nans
            pupilFitError = nan;
            pupilEllipseOnImagePlane=nan(1,5);
        else
            % For each position on the perimeter of the pupil, determine
            % which point (front or back) is farther from the center of the
            % ellipse. These will be the hidden points
            centerDistance = sqrt(sum(((...
                bsxfun(@minus,...
                mean([frontEllipse(1:2);backEllipse(1:2)]),...
                imagePoints(logical(pupilPerimIdxFront+pupilPerimIdxBack),:))...
                ).^2),2));
            hideBack = (centerDistance(1:nStopPerimPoints)-centerDistance(nStopPerimPoints+1:nStopPerimPoints*2))<0;
            % Identify the set of non-hidden pupil perimeter points
            backStopVisibleIdx = pupilPerimIdxBack;
            backStopVisibleIdx(pupilPerimIdxBack)=~hideBack;
            frontStopVisibleIdx = pupilPerimIdxFront;
            frontStopVisibleIdx(pupilPerimIdxFront)=hideBack;
            pupilPerimIdx = or(backStopVisibleIdx,frontStopVisibleIdx);
            % Remove those pupil perimeter points that have had poor ray
            % tracing
            pupilPerimIdx = and(pupilPerimIdx,targetIntersectError<p.Results.rayTraceErrorThreshold);
            % Fit the ellipse
            [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
            % Update the pointLabels to indicate the hidden points
            pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)),'_hidden');
            pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)),'_hidden');
        end
    end
else
    pupilFitError = nan;
    pupilEllipseOnImagePlane=nan(1,5);
end


end % projectModelEye