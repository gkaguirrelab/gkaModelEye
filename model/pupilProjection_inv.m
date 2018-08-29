function [eyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, varargin)
% Determine the eyePose corresponding to an observed entrance pupil ellipse
%
% Syntax:
%  [eyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry)
%
% Description:
%	Given the sceneGeometry and an ellipse on the image plane, this routine
%   finds parameters of the rotation of the eye and pupil radius that can
%   best account for the parameters of the ellipse.
%
% Notes:
%   Rotations - Eye rotations are given as azimuth and elevations in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds to
%   the position of the eye when a line that connects the center of
%   rotation of the eye with the center of the pupil is normal to the image
%   plane. Positive rotations correspond to rightward, upward, translation
%   of the pupil center in the image.
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
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
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
%                           torsion, pupil radius]. The default values here
%                           represent the physical limits of the projection
%                           model for azimuth, elevation, and pupil radius.
%                           Torsion is constrained to zero by default.
%  'centerErrorThresh'    - Scalar. The first of two stopping point
%                           criteria for the search.
%  'constraintTolerance'  - Scalar. The second of two stopping point
%                           criteria for the search. If passed, this value
%                           will over-ride the value in the sceneGeometry
%                           structure.
%  'repeatSearchThresh'   - Scalar. If a centerError output value is 
%                           obtained that is greater than this threshold,
%                           then a repeat search across eyePose values will
%                           be conducted to account for the possibility
%                           that the solution obtained was a local minimum.
%  'nMaxSearches'         - Scalar. The maximum number of searches that the
%                           routine will conduct as it attempts to avoid
%                           local minima.
%  'searchCount'          - Scalar. The number of searches that have been
%                           conducted so far. This value is modified as
%                           pupilProjection_inv calls itself recursively.
%                           It should not be modified by the external
%                           calling function.
%
% Outputs:
%   eyePose               - A 1x4 vector with values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius is in mm.
%   bestMatchEllipseOnImagePlane - A 1x5 vector that contains the
%                           parameters of pupil ellipse on the image plane
%                           cast in transparent form. This is the output of
%                           the pupilProjection_fwd model for the
%                           sceneGeometry and the found eyePose.
%   centerError           - Scalar. The Euclidean distance (in pixels)
%                           between the [x, y] center of the
%                           pupilEllipseOnImagePlane and the center of the
%                           bestMatchEllipseOnImagePlane.
%   shapeError            - Scalar. The proportion of error in fitting
%                           ellipse shape, range 0-1.
%   areaError             - Scalar. The proportion of error in fitting
%                           ellipse area; unbounded around zero.
%
% Examples:
%{
    %% Test if we can find the eyePose for image ellipse
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define in eyePoses the azimuth, elevation, torsion, and pupil radius
    eyePose = [10 10 0 2];
    % Obtain the pupil ellipse parameters in transparent format
    pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry);
    % Recover the eye pose from the ellipse
    inverseEyePose = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry);
    % Report the difference between the input and recovered eyePose
    fprintf('Test if the absolute error in the eye pose recovered by pupilProjection_inv is less than 1 percent\n');
    assert(max(abs(eyePose - inverseEyePose)./eyePose) < 1e-2)
%}
%{
    %% Calculate the time required for the inverse projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Generate ellipses for some randomly selected eye poses
    nPoses = 20;
    eyePoses=[(rand(nPoses,1)-0.5)*20, (rand(nPoses,1)-0.5)*10, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    for pp = 1:nPoses
    	ellipseParams(pp,:) = pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    fprintf('\nTime to compute inverse projection model (average over %d projections):\n',nPoses);
    tic
    for pp = 1:nPoses
    	pupilProjection_inv(ellipseParams(pp,:),sceneGeometry);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing pre-compiled ray tracing: %4.2f msecs.\n',msecPerModel);
%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('pupilEllipseOnImagePlane',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.5],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('centerErrorThresh',1e-4,@isnumeric);
p.addParameter('constraintTolerance',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('repeatSearchThresh',1.0,@isnumeric);
p.addParameter('nMaxSearches',3,@isnumeric);
p.addParameter('searchCount',1,@isnumeric);

% Parse and check the parameters
p.parse(pupilEllipseOnImagePlane, sceneGeometry, varargin{:});


%% Check inputs
% Handle an immediate exit
if isempty(pupilEllipseOnImagePlane)
    centerError=NaN;
    return
end

% Issue a warning if the bounds do not fully constrain at least one eye
% rotation parameter. This is because there are multiple combinations of
% the three axis rotations that can bring an eye to a destination.
% Typically, the torsion will be constrained with upper and lower bounds of
% zero, reflecting Listing*s Law.
if sum((p.Results.eyePoseUB(1:3) - p.Results.eyePoseLB(1:3))==0) < 1
    warning('pupilProjection_inv:underconstrainedSearch','The inverse search across possible eye rotations is underconstrained');
end

%% Assemble bounds and x0
% Because ellipses are symmetric about their axes, a given eccentricity and
% theta of an ellipse is consistent with two possible eye rotation
% solutions. We finesse this ambiguity by bounding the possible azimuth and
% elevation solutions depending upon the quadrant in which the center of
% the ellipse falls, relative to the center of projection derived from the
% scene geometry.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

% Identify the center of projection
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
CoP = [rotationCenterEllipse(1),rotationCenterEllipse(2)];

% Set the bounds on the eyePose based upon the quadrant of the ellipse
% center. We provide 5 degrees of wiggle in the fit around zero.
if pupilEllipseOnImagePlane(1) < CoP(1)
    eyePoseUB(1) = 5;
else
    eyePoseLB(1) = -5;
end
if pupilEllipseOnImagePlane(2) > CoP(2)
    eyePoseUB(2) = 5;
else
    eyePoseLB(2) = -5;
end

% If x0 is undefined, we make a guess based upon the location and size of
% the ellipse
if isempty(p.Results.x0)
    % Probe the forward model to determine how many pixels of change in the
    % location of the pupil ellipse correspond to one degree of rotation.
    % Omit ray-tracing to save time as it has minimal effect upon the
    % position of the center of the ellipse.
    probeEllipse=pupilProjection_fwd([1 0 0 2],sceneGeometry);
    pixelsPerDeg = probeEllipse(1)-CoP(1);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Torsion is set to
    % zero
    x0(1) = ((pupilEllipseOnImagePlane(1) - CoP(1))/pixelsPerDeg);
    x0(2) = ((CoP(2) - pupilEllipseOnImagePlane(2))/pixelsPerDeg);
    x0(3) = 0;
    
    % Force the angles of the x0 guess to be within bounds
    x0=min([eyePoseUB(1:3); x0]);
    x0=max([eyePoseLB(1:3); x0]);    
    
    % Estimate the pupil radius in pixels, accounting for the eccentricity
    % of the ellipse in the image plane
    ellipseAspectRatio = sqrt(1 - (pupilEllipseOnImagePlane(4)^2));
    pupilRadiusPixels = sqrt(pupilEllipseOnImagePlane(3) / (pi * ellipseAspectRatio));
    
    % Probe the forward model at the estimated pose angles to estimate the
    % pupil radius.
    probeEllipse=pupilProjection_fwd([x0(1) x0(2) x0(3) 2], sceneGeometry);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    x0(4) = pupilRadiusPixels/pixelsPerMM;
    
    % If the absolute value of an estimated angle is less than 2 degrees,
    % set the value to close to zero. This is done as fmincon seems to
    % avoid solutions exactly at zero, and this kludge fixes that behavior.
    x0(abs(x0(1:3))<2) = 1e-6;
    
    % Ensure that x0 lies within the bounds with a bit of headroom so that
    % the solver does not get stuck up against a bound.
    boundHeadroom = (eyePoseUB - eyePoseLB)*0.001;
    x0=min([eyePoseUB-boundHeadroom; x0]);
    x0=max([eyePoseLB+boundHeadroom; x0]);
else
    x0 = p.Results.x0;
end

% Generate a random x0 starting point if x0 contains an Inf flag
if any(isinf(x0))
    x0 = (eyePoseUB-eyePoseLB).*rand(1,4)+eyePoseLB;
end


%% Perform the search
% We use nested functions for the objective and constraint so that the
% forward pupil projection is only computed once for each iteration of the
% fmincon solver. We use a nested outfun that implements a custom stopping
% rule that speeds the search.

% Define variables used in the nested functions
targetEllipse = pupilEllipseOnImagePlane; % the target ellipse params
centerErrorThresh = p.Results.centerErrorThresh;
lastFVal = realmax;
bestFVal = realmax;
xLast = []; % Last place pupilProjection_fwd was called
xBest = []; % The x with the lowest objective function value that meets
            % the constraint tolerance
shapeErrorAtLast = 0;
shapeErrorAtBest = 0;
areaErrorAtLast = 0;
areaErrorAtBest = 0;
ellipseAtLast = []; % pupilProjection_fwd result at xLast
ellipseAtBest = []; % pupilProjection_fwd result at xBest

% Obtain the constraintTolerance
if isempty(p.Results.constraintTolerance)
    constraintTolerance = sceneGeometry.constraintTolerance;
else
    constraintTolerance = p.Results.constraintTolerance;
end

% Define search options. We use the sqp algorithm as it is content to find
% solutions equal to the upper or lower bounds.
options = optimoptions(@fmincon,...
    'Display','off', ...
    'Algorithm','sqp',...
    'OutputFcn',@outfun, ...
    'constraintTolerance',constraintTolerance);

% Call fmincon. Note the use of anonymous functions for the objective and
% constraint, which nested below.
fmincon(@objfun, x0, [], [], [], [], eyePoseLB, eyePoseUB, @constr, options);

    % Nested objective function
    function fval = objfun(x)
        % Check if computation is necessary
        if ~isequal(x,xLast)
            ellipseAtLast = pupilProjection_fwd(x, sceneGeometry);
            xLast = x;
        end
        % Compute objective function as Euclidean distance in the target
        % and candidate ellipse centers
        if any(isnan(ellipseAtLast))
            fval = realmax;
        else
            fval = sqrt((targetEllipse(1) - ellipseAtLast(1))^2 + ...
                (targetEllipse(2) - ellipseAtLast(2))^2);
        end
    end

    % Nested constraint function
    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) 
            ellipseAtLast = pupilProjection_fwd(x, sceneGeometry);
            xLast = x;
        end
        
        if any(isnan(ellipseAtLast))
            c = realmax;
            ceq = realmax;
            shapeErrorAtLast = c;
            areaErrorAtLast = ceq;
        else
            % c:
            % The theta and eccentricity of an ellipse can be described as
            % a point in polar coordinates. We express the constraint as
            % the vector distance between these points. Direct minimization
            % of differences in theta is a poor constraint, as differences
            % in theta have reduced meaning at small eccentricities.
            % Because ellipses are symmetric, theta spans the range of
            % 0:pi. Therefore, the theta value is doubled prior to
            % conversion to Cartesian coordinates so that the space wraps
            % at the 0 - pi transition point. Eccentricity has a value
            % ranging from zero (circular) to 1 (a fully flattened
            % ellipse). We linearize the eccentricity value so that the
            % error metric is sensitive to small differences in
            % eccentricity. The ceq value is divided by 2, so that the
            % largest possible error is unity.
            
            thetaT = targetEllipse(5)*2;
            thetaC = ellipseAtLast(5)*2;
            rhoT = 1-sqrt(1-targetEllipse(4)^2);
            rhoC = 1-sqrt(1-ellipseAtLast(4)^2);
            
            c = sqrt(rhoT^2 + rhoC^2 - 2*rhoT*rhoC*cos(thetaT-thetaC))/2;
            shapeErrorAtLast = c;
            
            % ceq:
            % Proportional difference in ellipse areas
            ceq = abs(targetEllipse(3) - ellipseAtLast(3))/targetEllipse(3);
            areaErrorAtLast = ceq;
        end
    end

    % Nested output function
    function stop = outfun(~,optimValues,state)
        stop = false;        
        switch state
            case 'init'
                lastFVal = optimValues.fval;
            case 'iter'
                lastFVal = optimValues.fval;
                % Store the current best value for x. This is done as we
                % observe that fmincon can move away from the best solution
                % when azimuth and elevation are close to zero. This
                % behavior has been seen by others:
                %   https://groups.google.com/forum/#!topic/comp.soft-sys.matlab/SuNzbhEun1Y
                if lastFVal < bestFVal
                    bestFVal = lastFVal;
                    xBest = xLast;
                    shapeErrorAtBest = shapeErrorAtLast;
                    areaErrorAtBest = areaErrorAtLast;
                    ellipseAtBest = ellipseAtLast;
                end
                % Test if we are done the search
                if lastFVal < centerErrorThresh && ...
                        optimValues.constrviolation < constraintTolerance
                    stop = true;
                end
            case 'done'
                % Unused
            otherwise
        end
    end


% Use the best solution seen by fmincon. This includes the eyePose, the
% parameters of the best fitting ellipse on the image plane, and the
% errors.
if isempty(xBest)
    eyePose = xLast;
    bestMatchEllipseOnImagePlane = ellipseAtLast;
    centerError = bestFVal;
    shapeError = shapeErrorAtLast;
    areaError = areaErrorAtLast;
else
    eyePose = xBest;
    bestMatchEllipseOnImagePlane = ellipseAtBest;
    centerError = bestFVal;
    shapeError = shapeErrorAtBest;
    areaError = areaErrorAtBest;
end

% If the solution has a centerError that is larger than repeatSearchThresh,
% we consider the possibility that the solution represents a local minimum.
% We repeat the search, passing a value close to the eyePose solution as
% x0. This process terminates when the search count exceeds nMaxSearches.
if centerError > p.Results.repeatSearchThresh && ...
        p.Results.searchCount < p.Results.nMaxSearches
    % Set the x0 for the next search to be the current eyePose solution
    x0 = eyePose;
    % Add a bit of offset to the azimuth and elevation of the x0 guess.
    % Scale this amount of offset by how many searches deep we are.
    x0(1:2) = x0(1:2)+[0.1 0.1]./p.Results.searchCount;
    % Execute a recursive call to pupilProjection_inv, supplying x0, and
    % the set of key-value pairs. Importantly, the searchCount value is
    % iterated.
    [eyePose_r, bestMatchEllipseOnImagePlane_r, centerError_r, shapeError_r, areaError_r] = ...
        pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, ...
        'x0',x0, ...
        'eyePoseLB',p.Results.eyePoseLB, ...
        'eyePoseUB',p.Results.eyePoseUB, ...
        'centerErrorThresh',p.Results.centerErrorThresh, ...
        'constraintTolerance',p.Results.constraintTolerance, ...
        'repeatSearchThresh', p.Results.repeatSearchThresh, ...
        'searchCount', p.Results.searchCount+1);
    % If the recursive search was better than the current search, save the
    % results. The best result across all recursive searches is thus
    % retained if the recusive stack collapses from exceeding the maximum
    % allowed number of searches.
    if centerError_r < centerError
        eyePose = eyePose_r;
        bestMatchEllipseOnImagePlane = bestMatchEllipseOnImagePlane_r;
        centerError = centerError_r;
        shapeError = shapeError_r;
        areaError = areaError_r;
    end
end

end % function -- pupilProjection_inv



