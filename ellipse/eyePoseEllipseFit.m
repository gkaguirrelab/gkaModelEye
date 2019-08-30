function [eyePose, RMSE, fittedEllipse, fitAtBound, nSearches] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Syntax:
%  [eyePose, RMSE, fittedEllipse, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry)
%
% Description:
%   The routine fits points on the image plane based upon the eye
%   parameters (azimuth, elevation, torsion, stop radius) that produce the
%   best fitting ellipse projected according to sceneGeometry.
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
%                           starting point will be estimated either from
%                           the polyModel field of the sceneGeometry or
%                           from the coordinates of the ellipse center.
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
    targetEllipse = pupilProjection_fwd(eyePose,sceneGeometry);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 5, 0 );
    [recoveredEyePose,RMSE,recoveredEllipse,fitAtBound,nSearches] = eyePoseEllipseFit(Xp, Yp, sceneGeometry);
    % Test that the recovered eye pose has no more than 0.01% error
    assert(max(abs((eyePose-recoveredEyePose)./eyePose)) < 1e-4)
%}


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('rmseThresh',1e-2,@isscalar);
p.addParameter('nMaxSearches',5,@isscalar);

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


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

%% Set bounds
% Identify the center of projection.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

% Clear the case of nans in the input
if any(isnan(eyePoseLB)) || any(isnan(eyePoseUB)) || any(any(isnan(p.Results.x0)))
    return
end


%% Define x0
% Check if x0 was passed
if isempty(p.Results.x0)
    % Define x0
    x0 = zeros(1,4);
    
    % Construct an x0 guess by probing the forward model. First identify
    % the center of projection
    rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
    CoP = [rotationCenterEllipse(1),rotationCenterEllipse(2)];
    
    meanXp = mean(Xp);
    meanYp = mean(Yp);
    
    % Probe the forward model to determine how many pixels of change in
    % the location of the pupil ellipse correspond to one degree of
    % rotation.
    probeEllipse=pupilProjection_fwd([1 1 0 2],sceneGeometry);
    pixelsPerDegHorz = probeEllipse(1)-CoP(1);
    pixelsPerDegVert = probeEllipse(2)-CoP(2);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection.
    x0(1) = ((meanXp - CoP(1))/pixelsPerDegHorz);
    x0(2) = ((meanYp - CoP(2))/pixelsPerDegVert);
    
    % Force the angles within bounds
    x0=[min([eyePoseUB(1:3); x0(1:3)]) 0];
    x0=[max([eyePoseLB(1:3); x0(1:3)]) 0];
    
    % Estimate the pupil radius in pixels
    pupilRadiusPixels = max([abs(max(Xp)-min(Xp)) abs(max(Yp)-min(Yp))])/2;
    
    % Probe the forward model at the estimated pose angles to estimate the
    % pupil radius.
    probeEllipse=pupilProjection_fwd([x0(1) x0(2) x0(3) 2], sceneGeometry);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    x0(4) = pupilRadiusPixels/pixelsPerMM;
    
    % Ensure that x0 lies within the bounds with a bit of headroom so that
    % the solver does not get stuck up against a bound.
    boundHeadroom = (eyePoseUB - eyePoseLB)*0.001;
    x0=min([eyePoseUB-boundHeadroom; x0]);
    x0=max([eyePoseLB+boundHeadroom; x0]);
else
    x0 = p.Results.x0;
end


%% Set up variables and functions for the search

% Define variables used in the nested functions
fittedEllipse = nan(1,5);

% define some search options
options = optimset('Display','off');

% Turn off warnings that can arise during the search
warningState = warning;
warning('off','pupilProjection_fwd:ellipseFitFailed');
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
    % Copy over the RMSE
    lastRMSE = RMSE;
    % Loop over the elements of the eyePose and search
    for ii = 1:length(eyePose)
        localObj = @(p) objfun(subP(ii,p,eyePose));
        [eyePose(ii), RMSE] = fminbnd(localObj, eyePoseLB(ii), eyePoseUB(ii), options);
    end
    % Check for termination conditions, which are any of:
    %   objective value is less than rmseThresh
    %   the change in objective val from last loop is less than rmseThresh
    %   we have used up all of our searches
    if RMSE < p.Results.rmseThresh || (lastRMSE-RMSE) < p.Results.rmseThresh || nSearches == p.Results.nMaxSearches
        searchingFlag = false;
    end
end % while
    function fVal = objfun(x)
        % Obtain the entrance pupil ellipse for this eyePose
        fittedEllipse = pupilProjection_fwd(x, sceneGeometry);
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
    end % local objective function


% Check if the fit is at a boundary for any parameter that is not locked
notLocked = eyePoseLB ~= eyePoseUB;
fitAtBound = any([any(abs(eyePose(notLocked)-eyePoseLB(notLocked))<1e-4) any(abs(eyePose(notLocked)-eyePoseUB(notLocked))<1e-4)]);

% Restore the warning state
warning(warningState);

end % eyePoseEllipseFit


