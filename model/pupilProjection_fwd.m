function [pupilEllipseOnImagePlane, imagePoints, worldPoints, eyePoints, pointLabels, nodalPointIntersectError, pupilFitError] = pupilProjection_fwd(eyePose, sceneGeometry, varargin)
% Obtain the parameters of the entrance pupil ellipse on the image plane
%
% Syntax:
%  [pupilEllipseOnImagePlane, imagePoints, worldPoints, eyePoints, pointLabels] = pupilProjection_fwd(eyePoses, sceneGeometry)
%
% Description:
%   Given the sceneGeometry this routine simulates an actual pupil on a
%   rotating eye and then measures the parameters of the ellipse (in
%   transparent format) of the entrance pupil in the image plane.
%
%   The forward model is a perspective projection of an anatomically
%   accurate eye, with points positioned behind the cornea subject to
%   refractive displacement. The projection incorporates the intrinsic
%   properties of the camera, including any radial lens distortion.
%
% Notes:
%   Rotations - Eye rotation is given as azimuth, elevation, and torsion in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds to
%   the position of the eye when a line that connects the center of
%   rotation of the eye with the center of the pupil is normal to the image
%   plane. Positive rotations correspond to rightward / upward (+x, -y)
%   translation of the pupil center in the image. Torsion of zero
%   corresponds to the torsion of the eye when it is in primary position.
%
%   Units - Eye rotations are in units of degrees. However, the units of
%   theta in the transparent ellipse parameters are radians. This is in
%   part to help us keep the two units separate conceptually.
%
% Inputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius in mm.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'fullEyeModelFlag'     - Logical. Determines if the full posterior and
%                           anterior chamber eye model will be created.
%  'nPupilPerimPoints'    - Scalar. The number of points that are 
%                           distributed around the pupil circle. A minimum
%                           of 5 is required to uniquely specify the image
%                           ellipse.
%  'pupilPerimPhase'      - Scalar. The phase (in radians) of the position
%                           of the pupil perimeter points.
%  'nIrisPerimPoints'     - The number of points that are distributed
%                           around the iris circle. A minimum of 5 is
%                           required to uniquely specify the image ellipse.
%  'posteriorChamberEllipsoidPoints' - The number of points that are on
%                           each latitude line of the posterior chamber
%                           ellipsoid. About 30 makes a nice image.
%  'anteriorChamberEllipsoidPoints' - The number of points that are on
%                           each longitude line of the anterior chamber
%                           ellipsoid. About 30 makes a nice image.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
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
%   eyePoints             - An nx3 matrix of the coordinates of the
%                           points of the eye model in the eyeWorld
%                           coordinate frame.
%   pointLabels           - An nx1 cell array that identifies each of the
%                           points, from the set {'pupilCenter',
%                           'irisCenter', 'aziRotationCenter',
%                           'eleRotationCenter', 'posteriorChamber',
%                           'irisPerimeter', 'pupilPerimeter',
%                           'anteriorChamber','cornealApex'}.
%   nodalPointIntersectError - A nx1 vector that contains the distance (in
%                           mm) between the nodal point of the camera and
%                           the intersection of a ray on the camera plane.
%                           The value is nan for points not subject to
%                           refraction by the cornea. All values will be
%                           nan sceneGeometry.refraction is empty.
%   pupilFitError         - The RMSE distance (in pixels) of the pupil
%                           perimeter points to the pupil ellipse.
%
% Examples:
%{
    %% Basic forward projection example
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define an eyePose with the azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Obtain the pupil ellipse parameters in transparent format
    pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry);
    % Test against cached result
    pupilEllipseOnImagePlaneCached = [0.027839132089926   0.022398638702744   1.547917599920448   0.000025771313824   0.000191221941167].*1e4;
    assert(max(abs(pupilEllipseOnImagePlane -  pupilEllipseOnImagePlaneCached)) < 1e-6)
%}
%{
    %% Display a 3D plot of a right eye
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Disable the refraction correction, as we are viewing the 3D model
    sceneGeometry.refraction = [];
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Perform the projection and request the full eye model
    [~, ~, worldPoints, ~, pointLabels] = pupilProjection_fwd(eyePose,sceneGeometry,'fullEyeModelFlag',true);
    % Define some settings for display
    eyePartLabels = {'aziRotationCenter', 'eleRotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeterFront' 'pupilPerimeterBack' 'anteriorChamber' 'cornealApex' 'fovea' 'opticDisc'};
    plotColors = {'>r' '^m' '.k' '*b' '*g' '*g' '.y' '*y' '*r' 'xk'};
    % Prepare a figure
    figure
    hold on
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)
    	idx = strcmp(pointLabels,eyePartLabels{pp});
        plot3(worldPoints(idx,1), worldPoints(idx,2), worldPoints(idx,3), plotColors{pp})
    end
    % Add the visual and optical axes
    opticalOrigin = worldPoints(strcmp(pointLabels,'opticalAxisOrigin'),:);
    foveaPoint = worldPoints(strcmp(pointLabels,'fovea'),:);
    nodalPoint = worldPoints(strcmp(pointLabels,'nodalPointRear'),:);
    visualAxis = [foveaPoint; nodalPoint];
    opticalAxis = [opticalOrigin; nodalPoint];
    plot3(visualAxis(:,1),visualAxis(:,2),visualAxis(:,3),'-r');
    plot3(opticalAxis(:,1),opticalAxis(:,2),opticalAxis(:,3),'-k');
    hold off
    axis equal
%}
%{
    %% Test the accuracy of the ellipse fit to the pupil boundary
    sceneGeometry=createSceneGeometry();
    aziVals = -60:5:60;
    pupilFitError = [];
    for aa = 1:length(aziVals)
        eyePose = [aziVals(aa) -3 0 3];
        [pupilEllipseOnImagePlane, ~, ~, ~, ~, ~, pupilFitError(aa)] = pupilProjection_fwd(eyePose, sceneGeometry,'nPupilPerimPoints',6);
    end
    figure
    plot(aziVals,pupilFitError,'.r');
%}
%{
    %% Calculate the time required for the forward projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Generate some randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*20, (rand(nPoses,1)-0.5)*10, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    fprintf('\nTime to compute forward projection model (average over %d projections):\n',nPoses);
    tic
    for pp = 1:nPoses
    	pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing pre-compiled ray tracing: %4.2f msecs.\n',msecPerModel);
%}



%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eyePose',@(x)(isnumeric(x) && all(size(x)==[1 4])));
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('fullEyeModelFlag',false,@islogical);
p.addParameter('nPupilPerimPoints',5,@(x)(isnumeric(x) && x>4));
p.addParameter('pupilPerimPhase',0,@isnumeric);
p.addParameter('nIrisPerimPoints',5,@isnumeric);
p.addParameter('posteriorChamberEllipsoidPoints',30,@isnumeric);
p.addParameter('anteriorChamberEllipsoidPoints',30,@isnumeric);

% parse
p.parse(eyePose, sceneGeometry, varargin{:})


%% Prepare variables
% Separate the eyePoses into individual variables
eyeAzimuth = eyePose(1);
eyeElevation = eyePose(2);
eyeTorsion = eyePose(3);
pupilRadius = eyePose(4);

% Store the number of pupil perimeter points
nPupilPerimPoints = p.Results.nPupilPerimPoints;


%% Define an eye in eye coordinates
% This coordinate frame is in mm units and has the dimensions (p1,p2,p3).
% The diagram is of a cartoon pupil, being viewed directly from the front.
%
% Coordinate [0,0,0] corresponds to the apex (front surface) of the cornea,
% with the model eye having the property of the optical and pupil axes of
% the eye being aligned. The first dimension is depth, and has a negative
% value towards the back of the eye.
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
% and positive values are more nasal. Positive values of p3 are downward,
% and negative values are upward


%% Define points around the elliptical actual pupil
% The eccentricity of the actual pupil is given by a stored function
actualPupilEccenFunc = str2func(sceneGeometry.eye.pupil.eccenFcnString);
% Determine the parameters of the ellipse that defines the actual pupil in
% the plane of the pupil. The absolute value of actualPupilEccenFunc gives
% the eccentricity. The theta of the theta of the actual pupil switches
% from horizontal to vertical when the pupil passes the circular radius
% point (0).
actualPupilEllipse = [sceneGeometry.eye.pupil.center(2) , ...
    sceneGeometry.eye.pupil.center(3), ...
    pi*pupilRadius^2, ...
    abs(actualPupilEccenFunc(pupilRadius)),...
    sceneGeometry.eye.pupil.thetas(1+(actualPupilEccenFunc(pupilRadius)>0))];
% Obtain the points on the perimeter of the ellipse
[p2p, p3p] = ellipsePerimeterPoints( actualPupilEllipse, nPupilPerimPoints, p.Results.pupilPerimPhase );
% Place these points into the eyeWorld coordinates, and create a front
% pupil perimeter
pupilFrontPoints(1:nPupilPerimPoints,3) = p3p;
pupilFrontPoints(1:nPupilPerimPoints,2) = p2p;
pupilFrontPoints(1:nPupilPerimPoints,1) = sceneGeometry.eye.pupil.center(1)+sceneGeometry.eye.iris.thickness/2;
eyePoints = pupilFrontPoints;
% Create labels for the pupilPerimeter points
tmpLabels = cell(nPupilPerimPoints, 1);
tmpLabels(:) = {'pupilPerimeterFront'};
pointLabels = tmpLabels;
% Add the back pupil perimeter
pupilBackPoints(1:nPupilPerimPoints,3) = p3p;
pupilBackPoints(1:nPupilPerimPoints,2) = p2p;
pupilBackPoints(1:nPupilPerimPoints,1) = sceneGeometry.eye.pupil.center(1)-sceneGeometry.eye.iris.thickness/2;
eyePoints = [eyePoints; pupilBackPoints];
tmpLabels(:) = {'pupilPerimeterBack'};
pointLabels = [pointLabels; tmpLabels];


%% Define full eye model
% If the fullEyeModel flag is set, then we will create a model of the
% posterior and anterior chambers of the eye.
if p.Results.fullEyeModelFlag
    
    % Add points for the pupil center, iris center, rotation centers,
    % origin of the optical axis, rear nodal point, and the fovea
    eyePoints = [eyePoints; sceneGeometry.eye.pupil.center];
    pointLabels = [pointLabels; 'pupilCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.iris.center];
    pointLabels = [pointLabels; 'irisCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.azi];
    pointLabels = [pointLabels; 'aziRotationCenter'];
    eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.ele];
    pointLabels = [pointLabels; 'eleRotationCenter'];
    eyePoints = [eyePoints; 0 0 0];
    pointLabels = [pointLabels; 'opticalAxisOrigin'];
    eyePoints = [eyePoints; sceneGeometry.eye.lens.nodalPoint];
    pointLabels = [pointLabels; 'nodalPoint'];
    eyePoints = [eyePoints; sceneGeometry.eye.posteriorChamber.fovea];
    pointLabels = [pointLabels; 'fovea'];
    eyePoints = [eyePoints; sceneGeometry.eye.posteriorChamber.opticDisc];
    pointLabels = [pointLabels; 'opticDisc'];
    
    % Define points around the perimeter of the iris
    nIrisPerimPoints = p.Results.nIrisPerimPoints;
    perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
    irisPoints(1:nIrisPerimPoints,3) = ...
        sin(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(3);
    irisPoints(1:nIrisPerimPoints,2) = ...
        cos(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(2);
    irisPoints(1:nIrisPerimPoints,1) = ...
        0 + sceneGeometry.eye.iris.center(1);
    
    % Add the points and labels
    eyePoints = [eyePoints; irisPoints];
    tmpLabels = cell(size(irisPoints,1), 1);
    tmpLabels(:) = {'irisPerimeter'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Create the anterior chamber ellipsoid
    [p1tmp, p2tmp, p3tmp] = ellipsoid( ...
        sceneGeometry.eye.cornea.front.center(1), ...
        sceneGeometry.eye.cornea.front.center(2), ...
        sceneGeometry.eye.cornea.front.center(3), ...
        sceneGeometry.eye.cornea.front.radii(1), ...
        sceneGeometry.eye.cornea.front.radii(2), ...
        sceneGeometry.eye.cornea.front.radii(3), ...
        p.Results.anteriorChamberEllipsoidPoints);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    anteriorChamberPoints=double(ansTmp.vertices);

    % Identify the index of the corneal apex
    [~,apexIdx]=max(anteriorChamberPoints(:,1));
    
    % Rotate the anteriorChamber points so that they reflect the difference
    % in the axis of the corneal ellipsoid w.r.t. the optical axis
    angles = sceneGeometry.eye.cornea.axis;
    R3 = [cosd(angles(1)) -sind(angles(1)) 0; sind(angles(1)) cosd(angles(1)) 0; 0 0 1];
    R2 = [cosd(-angles(2)) 0 sind(-angles(2)); 0 1 0; -sind(-angles(2)) 0 cosd(-angles(2))];
    R1 = [1 0 0; 0 cosd(angles(3)) -sind(angles(3)); 0 sind(angles(3)) cosd(angles(3))];
    anteriorChamberPoints = ((R1*R2*R3)*(anteriorChamberPoints-sceneGeometry.eye.cornea.front.center)')'+sceneGeometry.eye.cornea.front.center;
    
    % Save the corneal apex coordinates
    cornealApex = anteriorChamberPoints(apexIdx,:);
    
    % Retain those points that are anterior to the iris plane.
    retainIdx = anteriorChamberPoints(:,1) >= sceneGeometry.eye.iris.center(1);
    if all(~retainIdx)
        error('pupilProjection_fwd:pupilPlanePosition','The pupil plane is set in front of the corneal apea');
    end
    anteriorChamberPoints = anteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyePoints = [eyePoints; anteriorChamberPoints];
    tmpLabels = cell(size(anteriorChamberPoints,1), 1);
    tmpLabels(:) = {'anteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Add a entry for the corneal apex
    eyePoints = [eyePoints; cornealApex];
    pointLabels = [pointLabels; 'cornealApex'];
    
    % Create the posterior chamber ellipsoid. We switch dimensions here so
    % that the ellipsoid points have their poles at corneal apex and
    % posterior apex of the eye
    [p3tmp, p2tmp, p1tmp] = ellipsoid( ...
        sceneGeometry.eye.posteriorChamber.center(3), ...
        sceneGeometry.eye.posteriorChamber.center(2), ...
        sceneGeometry.eye.posteriorChamber.center(1), ...
        sceneGeometry.eye.posteriorChamber.radii(3), ...
        sceneGeometry.eye.posteriorChamber.radii(2), ...
        sceneGeometry.eye.posteriorChamber.radii(1), ...
        p.Results.posteriorChamberEllipsoidPoints);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    posteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are posterior to the iris plane, and have a
    % distance from the optical axis in the p2xp3 plane of greater than the
    % iris radius
    retainIdx = logical(...
        (posteriorChamberPoints(:,1) < sceneGeometry.eye.iris.center(1)) .* ...
        sqrt(posteriorChamberPoints(:,2).^2+posteriorChamberPoints(:,3).^2) > sceneGeometry.eye.iris.radius );
    if all(~retainIdx)
        error('pupilProjection_fwd:irisCenterPosition','The iris center is behind the center of the posterior chamber');
    end
    posteriorChamberPoints = posteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyePoints = [eyePoints; posteriorChamberPoints];
    tmpLabels = cell(size(posteriorChamberPoints,1), 1);
    tmpLabels(:) = {'posteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
end


%% Refract the eyeWorld points
% This steps accounts for the effect of corneal and corrective lens
% refraction upon the appearance of points from the eye.

% Define a variable to hold the calculated ray tracing errors
nodalPointIntersectError = nan(length(pointLabels),1);
% If we have a refraction field, proceed
if isfield(sceneGeometry,'refraction')
    % If this field is not set to empty, proceed    
    if ~isempty(sceneGeometry.refraction)
        % Assemble the static args for the virtualImageFunc
        args = {sceneGeometry.cameraPosition.translation, ...
                sceneGeometry.eye.rotationCenters, ...
                sceneGeometry.refraction.opticalSystem};
        % Identify the eyePoints subject to refraction by the cornea
        refractPointsIdx = find(...
            strcmp(pointLabels,'pupilPerimeterFront')+...
            strcmp(pointLabels,'pupilPerimeterBack')+...
            strcmp(pointLabels,'pupilCenter')+...
            strcmp(pointLabels,'irisPerimeter')+...
            strcmp(pointLabels,'irisCenter'));
        % Loop through the eyePoints that are to be refracted
        for ii=1:length(refractPointsIdx)
            % Get this eyeWorld point
            eyePoint=eyePoints(refractPointsIdx(ii),:);
            % Perform the computation using the passed function handle.
            [eyePoints(refractPointsIdx(ii),:), nodalPointIntersectError(refractPointsIdx(ii))] = ...
                sceneGeometry.refraction.handle(...
                eyePoint, eyePose, args{:});
        end
    end
end


%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to rightward, upward and clockwise (as seen
% from the perspective of the subject) eye movements
R.azi = [cosd(eyeAzimuth) -sind(eyeAzimuth) 0; sind(eyeAzimuth) cosd(eyeAzimuth) 0; 0 0 1];
R.ele = [cosd(-eyeElevation) 0 sind(-eyeElevation); 0 1 0; -sind(-eyeElevation) 0 cosd(-eyeElevation)];
R.tor = [1 0 0; 0 cosd(eyeTorsion) -sind(eyeTorsion); 0 sind(eyeTorsion) cosd(eyeTorsion)];


%% Apply the eye rotation
% This order (tor-ele-azi) corresponds to a head-fixed, extrinsic, rotation
% matrix. The reverse order (azi-ele-tor) would be an eye-fixed, intrinsic
% rotation matrix and would corresponds to the "Fick coordinate" scheme.
rotOrder = {'tor','ele','azi'};

% We shift the points to each rotation center, rotate, shift back, and
% repeat. Omit the eye rotation centers from this process and the center of
% projection. We must perform the rotation independently for each Euler
% angle to accomodate having rotation centers that differ by Euler angle.
rotatePointsIdx = ~contains(pointLabels,{'Rotation','centerOfProjection'});
for rr=1:3
    eyePoints(rotatePointsIdx,:) = ...
        (R.(rotOrder{rr})*(eyePoints(rotatePointsIdx,:)-sceneGeometry.eye.rotationCenters.(rotOrder{rr}))')'+sceneGeometry.eye.rotationCenters.(rotOrder{rr});
end

% If we are projecting a full eye model, label as hidden those posterior
% chamber points that are posterior to the most posterior of the centers of
% rotation of the eye, and thus would not be visible to the camera.
if p.Results.fullEyeModelFlag
    seenIdx = strcmp(pointLabels,'posteriorChamber') .* (eyePoints(:,1) >= min([sceneGeometry.eye.rotationCenters.azi(1) sceneGeometry.eye.rotationCenters.ele(1)]));
    seenIdx = logical(seenIdx + ~strcmp(pointLabels,'posteriorChamber'));
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
% pupil center when the line that connects the center of rotation of the
% eye and the pupil center are normal to the image plane.

% Re-arrange the eyePoints to transform to the world coordinate frame
worldPoints = eyePoints(:,[2 3 1]);


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
% for the position of the camera with respect to the world coordinate
% system. Only camera torsion is supported.
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
cameraExtrinsicMatrix = [transpose(cameraRotationMatrix) -transpose(cameraRotationMatrix)*sceneGeometry.cameraPosition.translation];

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
imagePointsNormalized = (imagePointsPreDistortion - [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]) ./ ...
    [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)];

% Distortion is proportional to distance from the center of projection on
% the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(1).*radialPosition.^2 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(2).*radialPosition.^4;

imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

% Place the distorted points back into the imagePoints vector
imagePoints = (imagePointsNormalizedDistorted .* [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]) +...
    [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)];


%% Obtain the pupil ellipse
% Proceed with fitting if we have a non-zero pupil radius
if eyePose(4) > 0
    % First obtain the ellipse fit to the front and back pupil perimeters
    pupilPerimIdxFront = strcmp(pointLabels,'pupilPerimeterFront');
    p1 = pupilEllipseFit(imagePoints(pupilPerimIdxFront,:));
    pupilPerimIdxBack = strcmp(pointLabels,'pupilPerimeterBack');
    p2 = pupilEllipseFit(imagePoints(pupilPerimIdxBack,:));
    if any(isnan(p1)) || any(isnan(p2))
        % If we are unable to fit an ellipse to the front or back pupil
        % perimeter, then exit with nans
        pupilFitError = nan;
        pupilEllipseOnImagePlane=nan(1,5);
    else
        % For each position on the perimeter of the pupil, determine which
        % point (front or back) is farther from the center of the ellipse
        % and then mark this point as hidden.
        centerDistance = sqrt(sum(((imagePoints(logical(pupilPerimIdxFront+pupilPerimIdxBack),:)-mean([p1(1:2);p2(1:2)])).^2),2));
        hideBack = (centerDistance(1:nPupilPerimPoints)-centerDistance(nPupilPerimPoints+1:nPupilPerimPoints*2))<0;
        idx = pupilPerimIdxBack;
        idx(pupilPerimIdxBack)=hideBack;
        pointLabels(idx) = strcat(pointLabels(idx),'_hidden');
        idx = pupilPerimIdxFront;
        idx(pupilPerimIdxFront)=~hideBack;
        pointLabels(idx) = strcat(pointLabels(idx),'_hidden');
        % Fit an ellipse to the non-hidden pupil perimeter points
        pupilPerimIdx = logical(strcmp(pointLabels,{'pupilPerimeterFront'}) + ...
            strcmp(pointLabels,{'pupilPerimeterBack'}));
        [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
    end
else
    pupilFitError = nan;
    pupilEllipseOnImagePlane=nan(1,5);
end


end % pupilProjection_fwd


%%%%% LOCAL FUNCTIONS

function [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints)

% Set up a variable to hold the pupil fit error
pupilFitError = nan;

% Before we try to fit the ellipse, make sure that the radius is not zero,
% and that there are at least 5 perimeter points that are non nan.
validPerimIdx = find(~any(isnan(imagePoints)')');

if ~isreal(imagePoints) || length(validPerimIdx)<5
    pupilEllipseOnImagePlane=nan(1,5);
else
    % Silence a warning that can arise regarding a nearly singular matrix
    warnState = warning;
    warning('off','MATLAB:nearlySingularMatrix');
    % We place the ellipse fit in a try-catch block, as the fit can fail
    % when the ellipse is so eccentric that it approaches a line
    try
        % Ellipse fitting with routine from the quadfit toolbox
        implicitEllipseParams = ellipsefit_direct( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2));
        % Convert the ellipse from implicit to transparent form
        pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
        % Place theta within the range of 0 to pi
        if pupilEllipseOnImagePlane(5) < 0
            pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
        end
        % Get the error of the ellipse fit to the pupil points
        pupilFitError = sqrt(nanmean(ellipsefit_distance( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2),ellipse_transparent2ex(pupilEllipseOnImagePlane)).^2));
    catch ME
        % If the ellipse fit direct fails because of an inability to fit an
        % ellipse to the provided points, return nans and issue a warning.
        pupilEllipseOnImagePlane=nan(1,5);
        pupilFitError = nan;
        switch ME.identifier
            case {'MATLAB:badsubscript','MATLAB:realsqrt:complexResult','MATLAB:expectedReal','MATLAB:quad2dproj:expectedFinite'}
                warning('pupilProjection_fwd:ellipseFitFailed','Could not fit a valid ellipse to the pupil points; returning nans.');
            otherwise
                warning('pupilProjection_fwd:ellipseFitUnknownError','Undefined error during ellipse fitting to pupil perimeter; returning nans.');
        end
    end % try-catch block
    warning(warnState);
end
end