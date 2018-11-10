function [eyePose, RMSE, ellipseParamsTransparent, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Syntax:
%  [eyePose, RMSE, ellipseParamsTransparent, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry)
%
% Description:
%   The routine fits points on the image plane based upon the eye
%   parameters (azimuth, elevation, torsion, pupil radius) that produce the
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
%  'x0'                   - Initial guess for the eyePose.
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%  'rmseThresh'           - Scalar that defines the stopping point for the
%                           search. The default value allows reconstruction
%                           of eyePose within 0.1% of the veridical,
%                           simulated value.
%
% Outputs:
%   eyePose               - A 1x4 matrix containing the best fitting eye
%                           parameters (azimuth, elevation, torsion, pupil
%                           radius)
%   RMSE                  - Root mean squared error of the distance of
%                           boundary point in the image to the fitted
%                           ellipse
%   ellipseParamsTransparent - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   fitAtBound            - Logical. Indicates if any of the returned
%                           ellipse parameters are at the upper or lower
%                           boundary.
%
% Examples:
%{
    %% Test if we can find the eyePose for a pupil perimeter ellipse
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define in eyePoses the azimuth, elevation, torsion, and pupil radius
    eyePose = [-25 25 0 2];
    % Obtain the pupil ellipse parameters in transparent format
    pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry);
    % Obtain boundary points for this ellipse. We need more than 5 boundary
    % points, as the pupil perimeter is not exactly elliptical
    [ Xp, Yp ] = ellipsePerimeterPoints( pupilEllipseOnImagePlane, 6 );
    % Recover the eye pose from the pupil boundary
    inverseEyePose = eyePoseEllipseFit(Xp, Yp, sceneGeometry);
    % Report the difference between the input and recovered eyePose
    fprintf('Test if there is less than 0.1 percent absolute error in the recovered eye pose.\n');
    fitError = abs((eyePose - inverseEyePose)./eyePose);
    fitError = max(fitError(~isnan(fitError)));
    assert(fitError < 1e-3);
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
    fprintf('\nTime to compute inverse projection model from pupil perimeter (average over %d projections):\n',nPoses);
    tic
    for pp = 1:nPoses
        [ Xp, Yp ] = ellipsePerimeterPoints( ellipseParams(pp,:), 6 );
        eyePoseEllipseFit(Xp, Yp, sceneGeometry);
    end
    msecPerModel = toc / nPoses * 1000;
    fprintf('\tUsing pre-compiled ray tracing: %4.2f msecs.\n',msecPerModel);
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
p.addParameter('rmseThresh',1e-2,@isnumeric);
p.addParameter('repeatSearchThresh',1.0,@isnumeric);
p.addParameter('searchCount',1,@isnumeric);
p.addParameter('nMaxSearches',3,@isnumeric);

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


% Set the return variables in the event of an error
eyePose = [nan nan nan nan];
RMSE = nan;
ellipseParamsTransparent = [nan nan nan nan nan];
fitAtBound = false;

% Issue a warning if the bounds do not fully constrain at least one eye
% rotation parameter. This is because there are multiple combinations of
% the three axis rotations that can bring an eye to a destination.
% Typically, the torsion will be constrained with upper and lower bounds of
% zero, reflecting Listing's Law.
if sum((p.Results.eyePoseUB(1:3) - p.Results.eyePoseLB(1:3))==0) < 1
    warning('eyePoseEllipseFit:underconstrainedSearch','The eye pose search across possible eye rotations is underconstrained');
end

%% Set bounds and x0
% Identify the center of projection.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

% Clear the case of nans in the input
if any(isnan(eyePoseLB)) || any(isnan(eyePoseUB)) || any(isnan(p.Results.x0)) 
    return
end

% Identify the center of projection
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
CoP = [rotationCenterEllipse(1),rotationCenterEllipse(2)];

meanXp = mean(Xp);
meanYp = mean(Yp);

% Set the bounds on the eyePose based upon the quadrant of the ellipse
% center. We provide a few degrees of wiggle in the fit around zero.
if meanXp < CoP(1)
    eyePoseUB(1) = 5;
else
    eyePoseLB(1) = -5;
end
if meanYp > CoP(2)
    eyePoseUB(2) = 5;
else
    eyePoseLB(2) = -5;
end

% If x0 is undefined, we make a guess based upon the location of the center
% of the points to be fit
if isempty(p.Results.x0)
    % Probe the forward model to determine how many pixels of change in the
    % location of the pupil ellipse correspond to one degree of rotation.
    probeEllipse=pupilProjection_fwd([1 0 0 2],sceneGeometry);
    pixelsPerDeg = probeEllipse(1)-CoP(1);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Torsion is set to
    % zero
    x0(1) = ((meanXp - CoP(1))/pixelsPerDeg);
    x0(2) = ((CoP(2) - meanYp)/pixelsPerDeg);
    x0(3) = 0;
    
    % Force the angles within bounds
    x0=min([eyePoseUB(1:3); x0]);
    x0=max([eyePoseLB(1:3); x0]);
    
    % Estimate the pupil radius in pixels
    pupilRadiusPixels = max([abs(max(Xp)-min(Xp)) abs(max(Yp)-min(Yp))])/2;
    
    % Probe the forward model at the estimated pose angles to
    % estimate the pupil radius.
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


% Define variables used in the nested functions
lastFVal = nan;
bestFVal = nan;
xLast = [nan nan nan nan]; % Last place pupilProjection_fwd was called
xBest = [nan nan nan nan]; % The x with the lowest objective function value
rmseThresh = p.Results.rmseThresh;
ellipseParamsTransparentLast = [nan nan nan nan nan];
ellipseParamsTransparentBest = [nan nan nan nan nan];

% define some search options
options = optimoptions(@fmincon,...
    'Display','off', ...
    'Algorithm','interior-point',...
    'OutputFcn',@outfun);

% Turn off warnings that can arise during the search
warningState = warning;
warning('off','pupilProjection_fwd:ellipseFitFailed');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% Clear the warning buffer
lastwarn('');

% Perform the search with nested objfun and outfun. We place the fmincon
% search within a try-catch block, as mysterious errors rarely occur deep
% within the fmincon routine, and I can't figure out why.
try
    fmincon(@objfun, x0, [], [], [], [], eyePoseLB, eyePoseUB, [], options);
catch
    eyePose = xBest;
    RMSE = bestFVal;
    ellipseParamsTransparent = ellipseParamsTransparentBest;
    % Check if the fit is at a boundary
    fitAtBound = any([any(eyePose==eyePoseLB) any(eyePose==eyePoseUB)]);
    return
end
    function fVal = objfun(x)
        xLast = x;
        % Obtain the entrance pupil ellipse for this eyePose
        ellipseParamsTransparentLast = pupilProjection_fwd(x, sceneGeometry);
        % Check for the case in which the transparentEllipse contains nan
        % values, which can arise if there were an insufficient number of
        % pupil border points remaining after refraction to define an
        % ellipse.
        if any(isnan(ellipseParamsTransparentLast))
            fVal = nan;
        else
            % This is the RMSE of the distance values of the boundary
            % points to the ellipse fit.
            explicitEllipse = ellipse_transparent2ex(ellipseParamsTransparentLast);
            if isempty(explicitEllipse)
                fVal = nan;
            else
                if any(isnan(explicitEllipse))
                    fVal = nan;
                else
                    fVal = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
                end
            end
        end
    end % local objective function
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
                if lastFVal < bestFVal || isnan(bestFVal)
                    bestFVal = lastFVal;
                    xBest = xLast;
                    ellipseParamsTransparentBest = ellipseParamsTransparentLast;
                end
                % Test if we are done the search
                if lastFVal < rmseThresh
                    stop = true;
                end
            case 'done'
                % Unused
            otherwise
        end
    end

% Store the best performing values
eyePose = xBest;
RMSE = bestFVal;
ellipseParamsTransparent = ellipseParamsTransparentBest;

% Check if the fit is at a boundary for any parameter that is not locked
notLocked = eyePoseLB ~= eyePoseUB;
fitAtBound = any([any(eyePose(notLocked)==eyePoseLB(notLocked)) any(eyePose(notLocked)==eyePoseUB(notLocked))]);

% Restore the warning state
warning(warningState);

% If a warning was received during the execution, exit here to avoid being
% stuck in a recursion loop of bad fits
[~,warnID] = lastwarn();
if ~isempty(warnID)
    % message = ['Received warning during eyePoseEllipseFit: ' warnID ' \n'];
    % fprintf(message);
    return
end

% If the solution has an RMSE that is larger than repeatSearchThresh, we
% consider the possibility that the solution represents a local minimum. We
% repeat the search, passing a value close to the eyePose solution as x0.
% This process terminates when the search count exceeds nMaxSearches.
if isnan(RMSE) || RMSE > p.Results.repeatSearchThresh
    if p.Results.searchCount < p.Results.nMaxSearches
        % Make sure that the search yielded an actual solution for the eyePose.
        % If not, simply re-use the x0 (with a small shift).
        if ~isempty(eyePose)
            x0 = eyePose;
        end
        x0(1:2) = x0(1:2)+[0.1 0.1]./p.Results.searchCount;
        [eyePose_r, RMSE_r, ellipseParamsTransparent_r, fitAtBound_r] = ...
            eyePoseEllipseFit(Xp, Yp, sceneGeometry, ...
            'x0',x0,...
            'eyePoseLB',p.Results.eyePoseLB,...
            'eyePoseUB',p.Results.eyePoseUB,...
            'repeatSearchThresh',p.Results.repeatSearchThresh, ...
            'searchCount', p.Results.searchCount+1);
        % Keep this solution if it is better
        if RMSE_r < RMSE
            eyePose = eyePose_r;
            RMSE = RMSE_r;
            ellipseParamsTransparent = ellipseParamsTransparent_r;
            fitAtBound = fitAtBound_r;
        end
    end
end

% If the eyePose variable is empty, then a fit could not be found. In this 
% situation return nans for the eyePose and nan for the RMSE
if isempty(eyePose)
    eyePose = [nan nan nan nan];
    RMSE = nan;
    ellipseParamsTransparent = [nan nan nan nan nan];
    fitAtBound = false;
end

end % eyePoseEllipseFit


