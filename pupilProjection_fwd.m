function [pupilEllipseOnImagePlane, imagePoints, sceneWorldPoints, eyeWorldPoints, pointLabels, nodalPointIntersectError, pupilFitError] = pupilProjection_fwd(eyePose, sceneGeometry, varargin)
% Project the pupil circle to an ellipse on the image plane
%
% Syntax:
%  [pupilEllipseOnImagePlane, imagePoints, sceneWorldPoints, eyeWorldPoints, pointLabels] = pupilProjection_fwd(eyePoses, sceneGeometry)
%
% Description:
%   Given the sceneGeometry this routine simulates an elliptical exit pupil
%   on a rotating eye and then measures the parameters of the ellipse (in
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
%   plane. Positive rotations correspond to rightward, upward, translation
%   of the pupil center in the image. Torsion of zero corresponds to the
%   torsion of the eye when it is in primary position.
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
%  'nPupilPerimPoints'    - The number of points that are distributed
%                           around the pupil circle. A minimum of 5 is
%                           required to uniquely specify the image ellipse.
%  'nIrisPerimPoints'     - The number of points that are distributed
%                           around the iris circle. A minimum of 5 is
%                           required to uniquely specify the image ellipse.
%  'posteriorChamberEllipsoidPoints' - The number of points that are on
%                           each latitude line of the posterior chamber
%                           ellipsoid. About 30 makes a nice image.
%  'anteriorChamberEllipsoidPoints' - The number of points that are on
%                           each longitude line of the anterior chamber
%                           ellipsoid. About 30 makes a nice image.
%  'removeOccultedPoints' - Logical. If set to true (default) the set of
%                           modeled points will be restricted to only those
%                           that are anterior to the center of rotation of
%                           the eye following rotation of the eye.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
%   imagePoints           - An nx2 matrix that specifies the x, y location
%                           on the image plane for each of the eyeWorld
%                           points.
%   sceneWorldPoints      - An nx3 matrix of the coordinates of the
%                           points of the eye model in the sceneWorld
%                           coordinate frame. If fullEyeModel is set to
%                           false, then there will only be points for the
%                           pupil perimeter. If fullEyeModel is true, then
%                           the entire model of ~1000 points will be
%                           returned.
%   eyeWorldPoints        - An nx3 matrix of the coordinates of the
%                           points of the eye model in the eyeWorld
%                           coordinate frame.
%   pointsLabels          - An nx1 cell array that identifies each of the
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
%                           nan sceneGeometry.virtualImageFunc is empty.
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
    % Test against 4/1/2018 cached result for eyePose [-10 5 0 3]
    assert(max(abs(pupilEllipseOnImagePlane -  [2.739864849789935e+02 2.215036804371633e+02 1.764122079488454e+04 0.193258910639646 2.151744146100150])) < 1e-6)
%}
%{
    %% Basic forward projection with a compiled virtualImageFunc
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
    % Define an eyePose with the azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Obtain the pupil ellipse parameters in transparent format
    pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry);
    % Test against 4/1/2018 cached result for eyePose [-10 5 0 3]
    assert(max(abs(pupilEllipseOnImagePlane -  [2.739864849789935e+02 2.215036804371633e+02 1.764122079488454e+04 0.193258910639646 2.151744146100150])) < 1e-6)
%}
%{
    %% Plot the pupil ellipse for various eye poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Prepare a figure
    figure
    sensorResolution=sceneGeometry.cameraIntrinsic.sensorResolution;
    blankFrame = zeros(fliplr(sensorResolution))+0.5;
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 sensorResolution(1)]);
    ylim([0 sensorResolution(2)]);
    % Loop over eye poses and plot
    for azi = -35:35:35
        for ele = -35:35:35
            eyePose = [azi ele 0 1];
            % Obtain the pupil ellipse parameters in transparent format
            pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry);
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilEllipseOnImagePlane));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            fimplicit(fh,[1, sensorResolution(1), 1, sensorResolution(2)],'Color', 'g','LineWidth',1);
            axis off;
        end
    end
    hold off
%}
%{
    %% Display a 2D image of a slightly myopic left eye wearing a contact
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry('eyeLaterality','left','sphericalAmetropia',-2,'contactLens',-2);
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 -5 0 3];
    % Perform the projection and request the full eye model
    [~, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd(eyePose,sceneGeometry,'fullEyeModelFlag',true);
    % Define some settings for display
    eyePartLabels = {'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber'};
    plotColors = {'.w' '.b' '*g' '.y'};
    sensorResolution=sceneGeometry.cameraIntrinsic.sensorResolution;
    blankFrame = zeros(fliplr(sensorResolution))+0.5;
    % Prepare a figure
    figure
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 sensorResolution(1)]);
    ylim([0 sensorResolution(2)]);
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)
    	idx = strcmp(pointLabels,eyePartLabels{pp});
        plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{pp})
    end
%}
%{
    %% Display a 3D plot of a right eye
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Disable the virtualImageFunc, as we are viewing the 3D model
    sceneGeometry.vitualImageFunc = [];
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Perform the projection and request the full eye model
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd(eyePose,sceneGeometry,'fullEyeModelFlag',true,'removeOccultedPoints',false);
    % Define some settings for display
    eyePartLabels = {'aziRotationCenter', 'eleRotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex'};
    plotColors = {'>r' '^m' '.k' '*b' '*g' '.y' '*y'};
    % Prepare a figure
    figure
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)
    	idx = strcmp(pointLabels,eyePartLabels{pp});
        plot3(sceneWorldPoints(idx,1), sceneWorldPoints(idx,2), sceneWorldPoints(idx,3), plotColors{pp})
        hold on
    end
    hold off
    axis equal
%}
%{
    %% Calculate the ray tracing error for some random poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Perform 100 forward projections with randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*30, (rand(nPoses,1)-0.5)*20, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    for pp = 1:nPoses
    	[~,~,~,~,~,nodalPointIntersectError(:,pp)]=pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    % Observe that the ray trace nodal error, while small, grows as a
    % function of the rotation of the eye.
    figure
    plot(sqrt(eyePoses(:,1).^2+eyePoses(:,2).^2),median(nodalPointIntersectError),'.r')
    xlabel('Euclidean rotation distance [deg]');
    ylabel('Ray trace nodal error [mm]');
%}
%{
    %% Calculate the time required for the forward projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Generate some randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*20, (rand(nPoses,1)-0.5)*10, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    clc
    fprintf('\nTime to compute forward projection model (average over %d projections):\n',nPoses);
    % Perform the forward projection without ray tracing
    sg = sceneGeometry;
    sg.virtualImageFunc = [];
    tic
    for pp = 1:nPoses
    	pupilProjection_fwd(eyePoses(pp,:),sg);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tWithout ray tracing: %4.2f msecs.\n',msecPerModel);
    % With ray tracing, using MATLAB routine
    sg = sceneGeometry;
    tic
    for pp = 1:nPoses
    	pupilProjection_fwd(eyePoses(pp,:),sg);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing MATLAB ray tracing: %4.2f msecs.\n',msecPerModel);
    % With ray tracing, using a compiled ray tracing routine
    sg = sceneGeometry;
    sg.virtualImageFunc = compileVirtualImageFunc(sg,'/tmp/demo_virtualImageFunc');
    tic
    for pp = 1:nPoses
    	pupilProjection_fwd(eyePoses(pp,:),sg);
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
p.addParameter('nIrisPerimPoints',5,@(x)(isnumeric(x) && x>4));
p.addParameter('posteriorChamberEllipsoidPoints',30,@isnumeric);
p.addParameter('anteriorChamberEllipsoidPoints',30,@isnumeric);
p.addParameter('removeOccultedPoints',true,@islogical);

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


%% Define an eye in eyeWorld coordinates
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
%  -  |        /     \
% p3  -  -----(   +   )-----
%  +  |        \_____/
%     v           |
%                 |
%
%           - <--p2--> +
%
% For the right eye, negative values on the p2 dimension are more temporal,
% and positive values are more nasal. Positive values of p3 are downward,
% and negative values are upward

%% Define points around the elliptical exit pupil
% The eccentricity of the exit pupil is given by a stored function
exitPupilEccenFunc = str2func(sceneGeometry.eye.pupil.eccenFcnString);
% Determine the parameters of the ellipse that defines the exit
% pupil in the plane of the pupil. The absolute value of exitPupilEccenFunc
% gives the eccentricity. The theta of the theta of the exit pupil switches
% from horizontal to vertical when the pupil passes the circular radius
% point (0).
exitPupilEllipse = [sceneGeometry.eye.pupil.center(2) , ...
    sceneGeometry.eye.pupil.center(3), ...
    pi*pupilRadius^2, ...
    abs(exitPupilEccenFunc(pupilRadius)),...
    sceneGeometry.eye.pupil.thetas(1+(exitPupilEccenFunc(pupilRadius)>0))];
% Obtain the points on the perimeter of the ellipse
[p2p, p3p] = ellipsePerimeterPoints( exitPupilEllipse, nPupilPerimPoints );
% Place these points into the eyeWorld coordinates
eyeWorldPoints(1:nPupilPerimPoints,3) = p3p;
eyeWorldPoints(1:nPupilPerimPoints,2) = p2p;
eyeWorldPoints(1:nPupilPerimPoints,1) = sceneGeometry.eye.pupil.center(1);
% Create labels for the pupilPerimeter points
tmpLabels = cell(nPupilPerimPoints, 1);
tmpLabels(:) = {'pupilPerimeter'};
pointLabels = tmpLabels;


%% Define full eye model
% If the fullEyeModel flag is set, then we will create a model of the
% posterior and anterior chambers of the eye.
if p.Results.fullEyeModelFlag
    
    % Add points for the center of the pupil, iris, and rotation
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.pupil.center];
    pointLabels = [pointLabels; 'pupilCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.iris.center];
    pointLabels = [pointLabels; 'irisCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.rotationCenters.azi];
    pointLabels = [pointLabels; 'aziRotationCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.rotationCenters.ele];
    pointLabels = [pointLabels; 'eleRotationCenter'];
    
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
    eyeWorldPoints = [eyeWorldPoints; irisPoints];
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
    anteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are anterior to the iris plane and not at a
    % greater radius in the p2xp3 plane than the iris.
    retainIdx = logical(...
        (anteriorChamberPoints(:,1) >= sceneGeometry.eye.iris.center(1)) .* ...
        (sqrt(anteriorChamberPoints(:,2).^2+anteriorChamberPoints(:,3).^2) <= sceneGeometry.eye.iris.radius) ...
        );
    if all(~retainIdx)
        error('pupilProjection_fwd:pupilPlanePosition','The pupil plane is set in front of the corneal apea');
    end
    anteriorChamberPoints = anteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; anteriorChamberPoints];
    tmpLabels = cell(size(anteriorChamberPoints,1), 1);
    tmpLabels(:) = {'anteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Add a point for the corneal apex
    cornealApex=[0 0 0];
    eyeWorldPoints = [eyeWorldPoints; cornealApex];
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
    eyeWorldPoints = [eyeWorldPoints; posteriorChamberPoints];
    tmpLabels = cell(size(posteriorChamberPoints,1), 1);
    tmpLabels(:) = {'posteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
end


%% Project the eyeWorld points to headWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (h1,h2,h3).
% The diagram is of a cartoon eye, being viewed directly from the front.
%
%  h1 values negative --> towards the head, positive towards the camera
%
%         h2
%    0,0 ---->
%     |
%  h3 |
%     v
%
%               |
%             __|__
%            /  _  \
%    -------(  (_)  )-------  h2 (horizontal axis of the head)
%            \_____/          rotation about h2 causes pure vertical
%               |             eye movement
%               |
%
%               h3
%   (vertical axis of the head)
%  rotation about h3 causes pure
%     horizontal eye movement
%
%
%
% Position [0,-,-] indicates the front surface of the eye.
% Position [-,0,0] indicates the h2 / h3 position of the pupil axis of
% the eye when it is normal to the image plane.
%
% We will convert from this coordinate frame to that of the camera scene
% later.


%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to leftward, downward and clockwise (as seen
% from the perspective of the subject) eye movements
R.azi = [cosd(eyeAzimuth) -sind(eyeAzimuth) 0; sind(eyeAzimuth) cosd(eyeAzimuth) 0; 0 0 1];
R.ele = [cosd(eyeElevation) 0 sind(eyeElevation); 0 1 0; -sind(eyeElevation) 0 cosd(eyeElevation)];
R.tor = [1 0 0; 0 cosd(eyeTorsion) -sind(eyeTorsion); 0 sind(eyeTorsion) cosd(eyeTorsion)];


%% Obtain the virtual image for the eyeWorld points
% This steps accounts for the effect of corneal and corrective lens
% refraction upon the appearance of points from the eye.

% Define a variable to hold the calculated ray tracing errors
nodalPointIntersectError = nan(length(pointLabels),1);
% If we have a virtualImageFunction field, proceed
if isfield(sceneGeometry,'virtualImageFunc')
    % If this field is not set to empty, proceed    
    if ~isempty(sceneGeometry.virtualImageFunc)
        % If the virtualImageFunc is a compiled routine (which we detect by
        % checking if the function handle evaluates with exist() to 3),
        % check that the that the optical system in the function is the
        % same as that in the passed sceneGeometry. For now, this check is
        % deactivated, as it adds 2 msecs per forward model calculation.
        %{
        if exist(func2str(sceneGeometry.virtualImageFunc.handle))==3
            if ~(sceneGeometry.opticalSystem==sceneGeometry.virtualImageFunc.opticalSystem)
                warning('pupilProjection_fwd:opticalSystemMismatch','The optical system used to build the virtual image function does not match that in the sceneGeometry');
            end
        end
        %}
        % Identify the eyeWorldPoints subject to refraction by the cornea
        refractPointsIdx = find(strcmp(pointLabels,'pupilPerimeter')+...
            strcmp(pointLabels,'pupilCenter')+...
            strcmp(pointLabels,'irisCenter'));
        % Loop through the eyeWorldPoints that are to be refracted
        for ii=1:length(refractPointsIdx)
            % Get this eyeWorld point
            eyeWorldPoint=eyeWorldPoints(refractPointsIdx(ii),:);
            % Perform the computation using the passed function handle.
            % This occurs within a try-catch block, as the point to be
            % refracted may experience total internal reflection. When this
            % happens, the routine exits with an error.
            try
                [eyeWorldPoints(refractPointsIdx(ii),:), nodalPointIntersectError(refractPointsIdx(ii))] = ...
                    sceneGeometry.virtualImageFunc.handle(...
                    eyeWorldPoint, eyePose, ...
                    sceneGeometry.virtualImageFunc.args{:});
            catch ME
                eyeWorldPoints(refractPointsIdx(ii),:) = nan;
                nodalPointIntersectError(refractPointsIdx(ii)) = inf;
                ME.identifier
                switch ME.identifier
                    case 'Coder:toolbox:ElFunDomainError'
                        warning('pupilProjection_fwd:rayTracingError','virtualImageFuncMex experienced an error, perhaps due to total internal reflection at a surface; returning nans for this point.');
                    otherwise
                        warning('pupilProjection_fwd:rayTracingUnknown',['Received the error ' ME.identifier ' during ray tracing; returning nans for this point.']);
                end
            end
        end
    end
end


%% Apply the eye rotation
% Copy the eyeWorld points into headWorld
headWorldPoints=eyeWorldPoints;

% This order (tor-ele-azi) corresponds to a head-fixed, extrinsic, rotation
% matrix. The reverse order (azi-ele-tor) would be an eye-fixed, intrinsic
% rotation matrix and would corresponds to the "Fick coordinate" scheme.
rotOrder = {'tor','ele','azi'};

% We shift the headWorld points to this rotation center, rotate, shift
% back, and repeat. Omit the eye rotation centers from this process.
rotatePointsIdx = ~contains(pointLabels,'Rotation');
for rr=1:3
    headWorldPoints(rotatePointsIdx,:) = ...
        (R.(rotOrder{rr})*(headWorldPoints(rotatePointsIdx,:)-sceneGeometry.eye.rotationCenters.(rotOrder{rr}))')'+sceneGeometry.eye.rotationCenters.(rotOrder{rr});
end

% If we are projecting a full eye model, and the 'removeOccultedPoints' is
% set to true, then remove those points that are posterior to the most
% posterior of the centers of rotation of the eye, and thus would not be
% visible to the camera.
if p.Results.fullEyeModelFlag && p.Results.removeOccultedPoints
    retainIdx = headWorldPoints(:,1) >= min([sceneGeometry.eye.rotationCenters.azi(1) sceneGeometry.eye.rotationCenters.ele(1)]);
    eyeWorldPoints = eyeWorldPoints(retainIdx,:);
    headWorldPoints = headWorldPoints(retainIdx,:);
    pointLabels = pointLabels(retainIdx);
end


%% Project the headWorld points to sceneWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (taken from Leszek Swirski), being
% viewed from above:
%
%   |
%   |    .-.
%   |   |   | <- Head
%   |   `^u^'
% Z |      :V <- Camera    (As seen from above)
%   |      :
%   |      :
%  \|/     o <- Target
%
%     ----------> X
%
% +X = right
% +Y = up
% +Z = front (towards the camera)
%
% The origin [0,0,0] corresponds to the front surface of the eye and the
% pupil center when the line that connects the center of rotation of the
% eye and the pupil center are normal to the image plane.

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
sceneWorldPoints = headWorldPoints(:,[2 3 1]);


%% Project the sceneWorld points to the image plane
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

% Create the extrinsicRotationMatrix. The model specifies only the camera
% rotation about the Z axis of the sceneWorld coordinate system.
extrinsicRotationMatrix = ...
    [cosd(sceneGeometry.cameraExtrinsic.rotationZ)	-sind(sceneGeometry.cameraExtrinsic.rotationZ)	0; ...
    sind(sceneGeometry.cameraExtrinsic.rotationZ)     cosd(sceneGeometry.cameraExtrinsic.rotationZ)     0; ...
    0                                       0                                       1];

% Create the projectionMatrix
projectionMatrix = ...
    sceneGeometry.cameraIntrinsic.matrix * ...
    [extrinsicRotationMatrix, ...
    sceneGeometry.cameraExtrinsic.translation];

% What is our total number of points to project?
nEyeWorldPoints = size(eyeWorldPoints,1);

% Project the sceneWorld points to the image plane and scale. The
% sceneWorld points have a column of ones added to support the
% multiplication with a combined rotation and translation matrix
tmpImagePoints=(projectionMatrix*[sceneWorldPoints, ones(nEyeWorldPoints,1)]')';
imagePointsPreDistortion=zeros(nEyeWorldPoints,2);
imagePointsPreDistortion(:,1) = ...
    tmpImagePoints(:,1)./tmpImagePoints(:,3);
imagePointsPreDistortion(:,2) = ...
    tmpImagePoints(:,2)./tmpImagePoints(:,3);


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the sensor optical center (aka
% principal point), and the coordinates are in world units. To apply this
% distortion to our image coordinate points, we subtract the optical
% center, and then divide by fx and fy from the intrinsic matrix.
imagePointsNormalized = (imagePointsPreDistortion - [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)]) ./ ...
    [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)];

% Distortion is proportional to distance from the center of the center of
% projection on the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(1).*radialPosition.^2 + ...
    sceneGeometry.cameraIntrinsic.radialDistortion(2).*radialPosition.^4;

imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

% Place the distorted points back into the imagePoints vector
imagePoints = (imagePointsNormalizedDistorted .* [sceneGeometry.cameraIntrinsic.matrix(1,1) sceneGeometry.cameraIntrinsic.matrix(2,2)]) +...
    [sceneGeometry.cameraIntrinsic.matrix(1,3) sceneGeometry.cameraIntrinsic.matrix(2,3)];


%% Fit an ellipse to the pupil points in the image plane
% Obtain the transparent ellipse params of the projection of the pupil
% circle on the image plane.
pupilPerimIdx = find(strcmp(pointLabels,'pupilPerimeter'));

% Set up a variable to hold the pupil fit error
pupilFitError = nan;

% Before we try to fit the ellipse, make sure that the radius is not zero,
% and that there are at least 5 perimeter points that are non nan.
validPerimIdx = find(~any(isnan(imagePoints(pupilPerimIdx,:))')');

if eyePose(4)==0 || ~isreal(imagePoints(pupilPerimIdx,:)) || length(validPerimIdx)<5
    pupilEllipseOnImagePlane=nan(1,5);
else
    % We place the ellipse fit in a try-catch block, as the fit can fail
    % when the ellipse is so eccentric that it approaches a line
    try
        % Ellipse fitting with routine from the quadfit toolbox
        implicitEllipseParams = ellipsefit_direct( imagePoints(pupilPerimIdx(validPerimIdx),1), imagePoints(pupilPerimIdx(validPerimIdx),2));
        % Convert the ellipse from implicit to transparent form
        pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
        % place theta within the range of 0 to pi
        if pupilEllipseOnImagePlane(5) < 0
            pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
        end
        % Get the error of the ellipse fit to the pupil points
        pupilFitError = sqrt(nanmean(ellipsefit_distance( imagePoints(pupilPerimIdx(validPerimIdx),1), imagePoints(pupilPerimIdx(validPerimIdx),2),ellipse_transparent2ex(pupilEllipseOnImagePlane)).^2));
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
    end
end

end % pupilProjection_fwd



