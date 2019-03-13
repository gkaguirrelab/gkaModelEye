function [eyePose, RMSE, fittedEllipse, fitAtBound] = pupilProjection_inv(targetEllipse, sceneGeometry, varargin)
% Determine the eyePose corresponding to an observed entrance pupil ellipse
%
% Syntax:
%  [eyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError] = pupilProjection_inv(targetEllipse, sceneGeometry)
%
% Description:
%	Given the sceneGeometry and an ellipse on the image plane, this routine
%   finds parameters of the rotation of the eye and stop radius that can
%   best account for the parameters of the ellipse. This is implemented by
%   defining a set of points on the perimeter of the ellipse and then
%   passing these points to be fit by the routine eyePoseEllipseFit.
%
% Notes:
%   Rotations - Eye rotations are given as azimuth and elevations in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds to
%   the position of the eye when the optical axis of the model eye
%   intersects the aperture stop of the observing camera. Positive
%   rotations correspond to rightward, upward, translation of the pupil
%   center in the image.
%
%   The default values set for the bounds on these rotation values reflect
%   the physical limits of the projection model. Tighter, biologically
%   informed constraints may be passed by the calling function. Note that
%   the search is underconstrained if there is freedom in the values to be
%   found for azimuth, elevation, and torsion. Indeed, Listing's Law
%   describes the tendency of the eye (for head-fixed saccades) to hold
%   torsion to zero when rotation the eye to a new location. Therefore, the
%   upper and lower bounds on torsion should generally be set to zero,
%   unless there is some specific desire to model this component under
%   different circumstances (e.g., peripheral nystagmus).
%
%   Units - Eye rotations are in units of degrees. However, the units of
%   theta in the transparent ellipse parameters are radians. This is in
%   part to help us keep the two units separate conceptually.
%
% Inputs:
%   targetEllipse         - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'x0'                   - Starting point of the search for the eyePose.
%                           If not defined, the starting point will be
%                           estimated from the coordinates of the ellipse
%                           center. If set to Inf, a random x0 will be
%                           selected within the eyePose bounds.
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePose [azimuth, elevation,
%                           torsion, stop radius]. The default values here
%                           represent the physical limits of the projection
%                           model for azimuth, elevation, and stop radius.
%                           Torsion is constrained to zero by default.
%  'rmseThresh'           - Scalar that defines the stopping point for the
%                           search. The default value allows reconstruction
%                           of eyePose within 0.1% of the veridical,
%                           simulated value.
%  'repeatSearchThresh'   - Scalar. If the RMSE output value obtained is
%                           greater than this threshold, then a repeat
%                           search across eyePose values will be conducted
%                           to account for the possibility that the
%                           solution obtained was a local minimum.
%  'nMaxSearches'         - Scalar. The maximum number of searches that the
%                           routine will conduct as it attempts to avoid
%                           local minima.
%
% Outputs:
%   eyePose               - A 1x4 vector with values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius is in mm.
%   RMSE                  - Root mean squared error of the distance of
%                           boundary point in the image to the fitted
%                           ellipse
%   fittedEllipse         - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   fitAtBound            - Logical. Indicates if any of the returned
%                           eyePose parameters are at the upper or lower
%                           boundary.
%
% Examples:
%{
    %% Test if we can find the eyePose for image ellipse
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define in eyePoses the azimuth, elevation, torsion, and pupil radius
    eyePose = [10 10 0 2];
    % Obtain the pupil ellipse parameters in transparent format
    targetEllipse = pupilProjection_fwd(eyePose,sceneGeometry);
    % Recover the eye pose from the ellipse
    inverseEyePose = pupilProjection_inv(targetEllipse, sceneGeometry);
    % Report the difference between the input and recovered eyePose
    fprintf('Test if the absolute error in the eye pose recovered by pupilProjection_inv is less than 1 percent\n');
    assert(max(abs(eyePose - inverseEyePose)./eyePose) < 1e-2)
%}
%{
    %% Calculate the time required for the inverse projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Generate ellipses for some randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*20, (rand(nPoses,1)-0.5)*10, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    for pp = 1:nPoses
    	ellipseParams(pp,:) = pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    fprintf('\nTime to compute inverse projection model (average over %d projections):\n',nPoses);
    tic
    recoveredEyePoses = []; RMSEvals = [];
    for pp = 1:nPoses
    	[recoveredEyePoses(pp,:),RMSEvals(pp)] = pupilProjection_inv(ellipseParams(pp,:),sceneGeometry);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing pre-compiled ray tracing: %4.2f msecs.\n',msecPerModel);
    fprintf('Max errors in azi, ele, torsion, and stop radius:\n');
    max(eyePoses-recoveredEyePoses)
    fprintf('median RMSE:\n');
    median(RMSEvals)
%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('targetEllipse',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[],@(x)(isnumeric(x) || isempty(x)));
p.addParameter('eyePoseUB',[],@(x)(isnumeric(x) || isempty(x)));
p.addParameter('rmseThresh',[],@(x)(isscalar(x) || isempty(x)));
p.addParameter('repeatSearchThresh',[],@(x)(isscalar(x) || isempty(x)));
p.addParameter('nMaxSearches',[],@(x)(isscalar(x) || isempty(x)));

% Parse and check the parameters
p.parse(targetEllipse, sceneGeometry, varargin{:});


%% Set the return variables in the event of an error
eyePose = [nan nan nan nan];
RMSE = nan;
fittedEllipse = [nan nan nan nan nan];
fitAtBound = false;


%% Check inputs
% Handle an immediate exit
if isempty(targetEllipse)
    return
end


%% Define perimeter points
[ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 7, 0 );


%% Perform the search
[eyePose, RMSE, fittedEllipse, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin{:});


end % function -- pupilProjection_inv



