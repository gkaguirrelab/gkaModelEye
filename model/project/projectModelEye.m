function [pupilEllipse, glintCoord, imagePoints, worldPoints, headPoints, eyePoints, pointLabels, targetIntersectError, pupilFitError] = projectModelEye(eyePose, sceneGeometry, options)
% Obtain the parameters of the entrance pupil ellipse on the image plane
%
% Syntax:
%  [pupilEllipse, glintCoord, imagePoints, worldPoints, headPoints, eyePoints, pointLabels, targetIntersectError, pupilFitError] = projectModelEye(eyePose, sceneGeometry, options)
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
%   degrees. These values correspond to degrees of rotation of the eye in
%   the Fick coordinate sysem. Azimuth, Elevation of [0,0] corresponds to
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
%   'cameraTrans'         - A 3x1 vector with values for [horizontal;
%                           vertical; depth] in units of mm. The values are
%                           relative to:
%                               sceneGeometry.cameraPosition.translation
%                           This key-value provides an easy way to update
%                           the camera translation.
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
%   pupilEllipse          - A 1x5 vector with the parameters of the
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
    pupilEllipse = projectModelEye(eyePose,sceneGeometry);
    % Test against cached result
    pupilEllipseCached = [ 0.027836556248511   0.022389506399564   1.545403286693142   0.000023300752734   0.000192318702774 ].*1e4;
    assert(max(abs(pupilEllipse -  pupilEllipseCached)) < 1e-4)
%}
%{
    %% Test the accuracy of the ellipse fit to the pupil boundary
    sceneGeometry=createSceneGeometry();
    aziVals = -60:5:60;
    pupilFitError = [];
    for aa = 1:length(aziVals)
        eyePose = [aziVals(aa) -3 0 3];
        [pupilEllipse, ~, ~, ~, ~, ~, ~, ~, pupilFitError(aa)] = projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16);
    end
    figure
    plot(aziVals,pupilFitError,'.r');
%}
%{
    %% Show the non-elliptical iris perimeter at extreme viewing angles
    sceneGeometry=createSceneGeometry();
    eyePose = [-65 0 0 3];
    [pupilEllipse, ~, imagePoints, ~, ~, ~, pointLabels] = ...
        projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16, ...
        'replaceReflectedPoints',true, ...
        'nIrisPerimPoints',16,'fullEyeModelFlag', true);
    figure
    idx = strcmp(pointLabels,'pupilPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'.r');
    hold on
    addTransparentEllipseToFigure(pupilEllipse);
    idx = strcmp(pointLabels,'irisPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'.b');
    [pupilEllipse, ~, imagePoints, ~, ~, ~, pointLabels] = ...
        projectModelEye(eyePose, sceneGeometry,'nStopPerimPoints',16, ...
        'replaceReflectedPoints',false, ...
        'nIrisPerimPoints',16,'fullEyeModelFlag', true);
    idx = strcmp(pointLabels,'pupilPerimeter');
    plot(imagePoints(idx,1),imagePoints(idx,2),'xr');
    addTransparentEllipseToFigure(pupilEllipse);
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

%% argument block
arguments
    eyePose (1, 4) double {mustBeNumeric}
    sceneGeometry struct
    options.cameraTrans (3, 1) double = [0; 0; 0]
    options.addPseudoTorsion logical = true
    options.fullEyeModelFlag logical = false
    options.nStopPerimPoints double {mustBeScalarOrEmpty, mustBeGreaterThan(options.nStopPerimPoints,4)} = 6
    options.stopPerimPhase double {mustBeScalarOrEmpty} = 0
    options.replaceReflectedPoints logical = false
    options.borderSearchPrecision double {mustBeScalarOrEmpty} = 0.01
    options.rayTraceErrorThreshold double {mustBeScalarOrEmpty} = 0.01
    options.nIrisPerimPoints double {mustBeScalarOrEmpty} = 5
    options.corneaMeshDensity double {mustBeScalarOrEmpty} = 23
    options.retinaMeshDensity double {mustBeScalarOrEmpty} = 30
    options.pupilRayFunc {isFuncHandleOrEmpty(options.pupilRayFunc)} = @findPupilRayMex
    options.glintRayFunc {isFuncHandleOrEmpty(options.glintRayFunc)} = @findGlintRayMex
end


%% Define some variables
% Set these in case we exit early from the routine
imagePoints = [];
pupilEllipse=nan(1,5);
pupilFitError = nan;
glintCoord = [];


%% Update camera translation
sceneGeometry.cameraPosition.translation = ...
        sceneGeometry.cameraPosition.translation + options.cameraTrans;


%% Apply pseudoTorsion
eyePose = ...
    addPseudoTorsion(sceneGeometry,eyePose,options);


%% Define the aperture stop perimeter
[eyePoints, pointLabels] = ...
    addStopPerimeter(sceneGeometry,eyePose,options);


%% Add the rest of the eye
[eyePoints, pointLabels] = ...
    addFullEyeModel(eyePoints,pointLabels,sceneGeometry,options);


%% Subject the eye points to refraction
[eyePoints, pointLabels, targetIntersectError] = ...
    refractEyePoints(eyePoints,pointLabels,sceneGeometry,eyePose,options);


%% Add glint(s)
[eyePoints, pointLabels] = ...
    addGlint(eyePoints,pointLabels,sceneGeometry,eyePose,options);


%% Apply eye rotation
[headPoints, pointLabels] = ...
    applyEyeRotation(eyePoints,pointLabels,sceneGeometry,eyePose,options);


%% Switch to world coordinates.
worldPoints = ...
    switchCoordinates(headPoints);


%% Check if we can proceed with projection
if ~isfield(sceneGeometry,'cameraPosition') || ~isfield(sceneGeometry,'cameraIntrinsic')
    return
end


%% Project the world coordinate points to the image plane
imagePointsPreDistortion = ...
    projectToImagePlane(worldPoints,sceneGeometry);


%% Apply radial lens distortion
imagePoints = ...
    applyRadialLensDistortion(imagePointsPreDistortion,sceneGeometry);


%% Obtain the pupilEllipse
[pupilEllipse,pupilFitError,pointLabels] = ...
    obtainImagePlaneEllipse(imagePoints,pointLabels,sceneGeometry,eyePose,options);


%% Obtain the glintCoord
glintCoord = ...
    obtainGlintCoord(imagePoints,pointLabels);


end % projectModelEye