function [eyePose, RMSE, fittedEllipse, fitAtBound, nSearches] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Syntax:
%  [eyePose, RMSE, fittedEllipse, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry)
%
% Description:
%   The routine fits the pupil perimeter points on the image plane based
%   upon the eye parameters (azimuth, elevation, torsion, stop radius) that
%   produce the best fitting ellipse projected according to sceneGeometry.
%
%   The search is constrained by the upper and lower bounds of the eyePose.
%   The default values specified here represent the physical boundaries of
%   the rotation model. Tighter, biologically informed constraints may be
%   passed by the calling function.
%
% Inputs:
%   Xp, Yp                - Vector of points to be fit
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'x0'                   - A 1x4 vector that provides starting points for
%                           the search for the eyePose. If not defined, the
%                           starting point will be estimated.
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePose [azimuth, elevation,
%                           torsion, stop radius]. The default values here
%                           represent the physical limits of the projection
%                           model for azimuth, elevation, and stop radius.
%                           Torsion is constrained to zero by default.
%  'rmseThresh'           - Scalar. The eyePose search will stop when the
%                           pupil perimeter points are fit with less than
%                           this error (in units of pixels). The search may
%                           also terminate with higher error if other
%                           stopping criteria are met.
%  'rmseThreshIter'       - Scalar. The eyePose search will stop when 
%                           sequential searches are producing a change in
%                           the RMSE of the fit than less than this value.
%  'eyePoseTol'           - Scalar. The eyePose values will be searched to
%                           within this level of precision.
%  'nMaxSearches'         - Scalar. The maximum number of searches that the
%                           routine will conduct as it attempts to avoid
%                           local minima.
%
% Outputs:
%   eyePose               - A 1x4 vector with values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius is in mm.
%   RMSE                  - Root mean squared error of the distance of
%                           boundary point in the image to the fitted
%                           ellipse
%   fittedEllipse         - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   fitAtBound            - Logical. Indicates if any of the returned
%                           eyePose parameters are at the upper or lower
%                           boundary.
%   nSearches             - The number of searches conducted.
%
% Examples:
%{
    sceneGeometry=createSceneGeometry();
    eyePose = [-15 10 0 2.5];
    [targetEllipse, glintCoord, imagePoints, ~, ~, ~, pointLabels] = projectModelEye(eyePose,sceneGeometry,'calcGlint',true);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 5, 0 );
    [recoveredEyePose,RMSE,recoveredEllipse,fitAtBound,nSearches] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, 'glintCoord', glintCoord);
    % Test that the recovered eye pose has no more than 0.01% error
    assert(max(abs((eyePose - recoveredEyePose)./eyePose)) < 1e-4)
%}


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('glintCoord',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('glintTol',1,@isscalar);
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('rmseThresh',[],@isscalar);
p.addParameter('rmseThreshIter',1e-2,@isscalar);
p.addParameter('eyePoseTol',1e-4,@isscalar);
p.addParameter('nMaxSearches',5,@isscalar);

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


%% Check inputs
% Initialize the return variables
eyePose = [nan nan nan nan];
RMSE = nan;
fittedEllipse = [nan nan nan nan nan];
fitAtBound = false;
nSearches = nan;

% Issue a warning if the bounds do not fully constrain at least one eye
% rotation parameter. This is because there are multiple combinations of
% the three axis rotations that can bring an eye to a destination.
% Typically, the torsion will be constrained with upper and lower bounds of
% zero, reflecting Listing's Law.
if sum((p.Results.eyePoseUB(1:3) - p.Results.eyePoseLB(1:3))==0) < 1
    warning('eyePoseEllipseFit:underconstrainedSearch','The eye pose search across possible eye rotations is underconstrained');
end

% If we have a glintCoord, make sure it is the right vector orientation
glintTol = p.Results.glintTol;
glintCoord = p.Results.glintCoord;
if ~isempty(glintCoord)
    if ~isequal(size(glintCoord),[1 2])
        error('eyePoseEllipseFit:glintCoordFormat','The optional glintCoord must be a 1x2 vector');
    end
end


%% Set bounds
% Identify the center of projection.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

% Identify the eyePose values that are free to vary in the search
notLocked = eyePoseLB ~= eyePoseUB;

% Clear the case of nans in the input
if any(isnan(eyePoseLB)) || any(isnan(eyePoseUB)) || any(any(isnan(p.Results.x0)))
    return
end


%% Obtain the unconstrained ellipse fit to the perimeter
[unconstrainedEllipse, unconstrainedRMSE] = pupilEllipseFit([Xp;Yp]);

% Various degenerate sets of
% perimeter points can cause the ellipse fit to fail and return nans
if any(isnan(unconstrainedEllipse))
    % The fit failed. Sythesize an ellipse vector that has as its center
    % the mean of the X and Y positions of the perimeter points.
    unconstrainedEllipse = nan(1,5);
    unconstrainedEllipse(1:2) = [mean(Xp) mean(Yp)];
    % Set the unconstrainedRMSE to an arbitrary, relatively large number.
    unconstrainedRMSE = 1.0;
end


%% Determine the search termination parameters
% Set the rmseThresh
if isempty(p.Results.rmseThresh)
    % The best that this routine could do would be to match the fit
    % provided by an unconstrained ellipse fit.
    rmseThresh = unconstrainedRMSE;
else
    rmseThresh = p.Results.rmseThresh;
end

% When the search is not improving the RMSE by more than this, it is time
% to stop.
rmseThreshIter = p.Results.rmseThreshIter;    

% Set the eyePose search tolerance.
eyePoseTol = p.Results.eyePoseTol;


%% Define x0
% Check if x0 was passed
if isempty(p.Results.x0)
    % Define x0
    x0 = zeros(1,4);
    
    % Construct an x0 guess by probing the forward model. First identify
    % the center of projection
    rotationCenterEllipse = projectModelEye([0 0 0 2], sceneGeometry);
    CoP = [rotationCenterEllipse(1),rotationCenterEllipse(2)];
    
    % Now the number of pixels that the pupil center is displaced from the
    % CoP. Handle the direction of eye rotation
    displacePix = (unconstrainedEllipse(1:2)-CoP);
    displaceScaled = displacePix ./ max(displacePix);
    displaceScaled(2) = -displaceScaled(2);
    
    % Probe the forward model to determine how many pixels of change in the
    % location of the pupil ellipse correspond to one degree of rotation.
    probeEllipse = projectModelEye([displaceScaled 0 2],sceneGeometry);
    pixelsPerDeg = abs(probeEllipse(1:2)-CoP)./abs(displaceScaled);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Need to make the
    % second value negative to match the direction of rotation convention
    % for elevation.
    x0(1) = displacePix(1)/pixelsPerDeg(1);
    x0(2) = -displacePix(2)/pixelsPerDeg(2);
    
    % Force the angles within bounds
    x0=[min([eyePoseUB(1:3); x0(1:3)]) 0];
    x0=[max([eyePoseLB(1:3); x0(1:3)]) 0];
    
    % Estimate the pupil radius in pixels
    pupilRadiusPixels = max([abs(max(Xp)-min(Xp)) abs(max(Yp)-min(Yp))])/2;
    
    % Probe the forward model at the estimated pose angles to estimate the
    % pupil radius.
    probeEllipse=projectModelEye([x0(1) x0(2) x0(3) 2], sceneGeometry);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    x0(4) = pupilRadiusPixels/pixelsPerMM;
    
    % Ensure that x0 lies within the bounds with a bit of headroom so that
    % the solver does not get stuck up against a bound.
    boundHeadroom = (eyePoseUB - eyePoseLB)*0.001;
    x0=min([eyePoseUB-boundHeadroom; x0]);
    x0=max([eyePoseLB+boundHeadroom; x0]);
    
    % Force any locked parameter to have the locked value
    x0(~notLocked) = eyePoseUB(~notLocked);
    
else
    x0 = p.Results.x0;
end


%% Set up variables and functions for the search

% Define variables used in the nested functions
fittedEllipse = nan(1,5);

% define some search options for fminbnd
options = optimset('fminbnd');
options.Display = 'off';
options.TolX = eyePoseTol; % eyePose search precision
options.MaxFunEvals = 50; % I think this value is a good speed-accuracy trade off, but have not carefully tested

% Turn off warnings that can arise during the search
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% Clear the warning buffer
lastwarn('');

% Set some parameters for the search
nSearches = 0;
searchingFlag = true;
eyePose = x0;
RMSE = realmax;

% Anonymous function to update one element of the eyePose at a time
subP = @(ii,p,x) [x(1:ii-1) p x(ii+1:end)];

% Enter a while loop that iteratively refines the eyePose until criteria
% are met. The objfun is nested.
while searchingFlag
    % Iterate the search count
    nSearches = nSearches+1;
    % Copy over the RMSE and eyePose
    lastRMSE = RMSE;
    lastEyePose = eyePose;
    % Loop over the elements of the eyePose and search
    for ii = 1:length(eyePose)
        % Only search over the eyePose parameters that are not fully
        % constrained by the bounds
        if notLocked(ii)
            localObj = @(p) objfun(subP(ii,p,eyePose));
            [eyePose(ii), RMSE] = fminbnd(localObj, eyePoseLB(ii), eyePoseUB(ii), options);
        end
    end
    % Check for termination conditions, which are any of:
    %   - objective value is within rmseThreshIter of rmseThresh
    %   - the change in objective val from last loop is less than rmseThreshIter
    %   - the change in eyePose values are all below eyePoseTol
    %   - we have used up all of our searches
    if ...
            abs(RMSE-rmseThresh) < rmseThreshIter || ...
            abs(lastRMSE-RMSE) < rmseThreshIter || ...
            all(abs(lastEyePose-eyePose) < eyePoseTol) || ...
            nSearches == p.Results.nMaxSearches
        searchingFlag = false;
    end
end % while
    function fVal = objfun(x)
        % Obtain the entrance pupil ellipse for this eyePose
        [fittedEllipse, fittedGlint, imagePoints, ~, ~, ~, pointLabels] = projectModelEye(x, sceneGeometry, 'calcGlint',true);
        % Check for the case in which the transparentEllipse contains nan
        % values, which can arise if there were an insufficient number of
        % pupil border points remaining after refraction to define an
        % ellipse.
        if any(isnan(fittedEllipse))
            % Set fVal to something arbitrarily large
            fVal = 1e6;
        else
            % This is the RMSE of the distance values of the boundary
            % points to the ellipse fit.
            explicitEllipse = ellipse_transparent2ex(fittedEllipse);
            if isempty(explicitEllipse)
                fVal = 1e6;
            else
                if any(isnan(explicitEllipse))
                    fVal = 1e6;
                else
                    fVal = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
                end
            end
        end
        % Check the match to the glint
        if ~isempty(glintCoord)
            if isempty(fittedGlint)
                fVal = 1e6;
            else
                fVal = fVal * (1+norm(glintCoord - fittedGlint));
            end
        end
    end % local objective function


% Check if the fit is at a boundary for any parameter that is not locked
fitAtBound = any([any(abs(eyePose(notLocked)-eyePoseLB(notLocked)) < eyePoseTol) any(abs(eyePose(notLocked)-eyePoseUB(notLocked)) < eyePoseTol)]);

% Restore the warning state
warning(warningState);

end % eyePoseEllipseFit


