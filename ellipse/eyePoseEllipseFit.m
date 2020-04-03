function [eyePose, RMSE, fittedEllipse, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin)
% Return the eye pose that best fits an ellipse to the pupil perimeter
%
% Syntax:
%  [eyePose, RMSE, fittedEllipse, fitAtBound] = eyePoseEllipseFit(Xp, Yp, sceneGeometry)
%
% Description:
%   The routine fits the pupil perimeter points on the image plane based
%   upon the eye parameters (azimuth, elevation, torsion, stop radius) that
%   produce the best fitting ellipse based upon the sceneGeometry.
%
%   If one or more glint coordinates are supplied, the eyePose used to fit
%   the pupil perimeter is constrained to also produce a modeled location
%   of the glint(s) that matches the supplied coordinates.
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
%   glintCoord            - A nx2 vector with the image coordinates of the
%                           n glints.
%  'x0'                   - A 1x4 vector that provides starting points for
%                           the search for the eyePose. If not defined, the
%                           starting point will be estimated.
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePose [azimuth, elevation,
%                           torsion, stop radius]. The default values here
%                           represent the physical limits of the projection
%                           model for azimuth, elevation, and stop radius.
%                           Torsion is constrained to zero by default.
%  'eyePoseEllipseFitFunEvals' - Scalar. The maximum number of evals for
%                           the fminsearch operation. The example below
%                           examines the trade-off between execution time
%                           and function accuracy for this parameter. I
%                           find that 50 evals takes ~1 second per eyePose
%                           calculation on a laptop, and produces excellent
%                           accuracy with respect to the imprecision in
%                           empirical data.
%  'eyePoseTol'           - Scalar. The eyePose values will be searched to
%                           within this level of precision. Doesn't have
%                           too much of an effect upon execution time.
%
% Outputs:
%   eyePose               - A 1x4 vector with values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           degrees, and stop radius is in mm.
%   RMSE                  - Root mean squared error of the distance of
%                           boundary points in the image to the fitted
%                           ellipse
%   fittedEllipse         - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   fitAtBound            - Logical. Indicates if any of the returned
%                           eyePose parameters are at the upper or lower
%                           boundary.
%
% Examples:
%{
    % Explore the trade-off between accuracy and run time
    sceneGeometry=createSceneGeometry();
    executionTime = []; errors = [];
    poseCount = 0;
    evalNum = [5 10 20 50 100 200];
    for azi = [-10 -.5 0.5 10]; for ele = [-10 -.5 0.5 10]
        eyePose = [azi ele 0 2.5];
        poseCount = poseCount+1;
        [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry);
        [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 6, 0 );
        for ee = 1:length(evalNum)
            tic
            recoveredEyePose = eyePoseEllipseFit(Xp, Yp, sceneGeometry,'glintCoord', glintCoord, 'eyePoseEllipseFitFunEvals', evalNum(ee));
            executionTime(poseCount,ee) = toc();
            % Measure the absolute error in recovering eye pose values
            errors(poseCount,ee,:) = abs(eyePose - recoveredEyePose);
        end
    end; end
    timeByNumEvals = mean(executionTime);
    errorByNumEvals = max(squeeze(mean(errors(:,:,1:2)))');
    figHandle = figure;
    left_color = [1 0 0]; right_color = [0 0 0];
    set(figHandle,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis right
    semilogx(evalNum,timeByNumEvals,'-k')
    ylabel('Execution time per pose [secs]');
    yyaxis left
    semilogx(evalNum,errorByNumEvals,'-r')
    ylabel('Eye pose error [deg]');
    xlabel('Number of search evaluations');
    title('Performance of eyePoseEllipseFit across max fun evals');
%}
%{
    % Examine the behavior of the routine at the eyePose bounds
    sceneGeometry=createSceneGeometry();
    eyePose = [40 10 0 2.5];
    [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 6, 0 );
    recoveredEyePose = eyePoseEllipseFit(Xp, Yp, sceneGeometry,'glintCoord', glintCoord, 'eyePoseUB', [30 30 0 4]);
%}


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addParameter('glintCoord',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('eyePoseEllipseFitFunEvals',50,@isscalar);
p.addParameter('eyePoseTol',1e-3,@isscalar);


% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


%% Check inputs
% Initialize the return variables
eyePose = [nan nan nan nan];
RMSE = nan;
fittedEllipse = [nan nan nan nan nan];
fitAtBound = false;

% Issue a warning if the bounds do not fully constrain at least one eye
% rotation parameter. This is because there are multiple combinations of
% the three axis rotations that can bring an eye to a destination.
% Typically, the torsion will be constrained with upper and lower bounds of
% zero, reflecting Listing's Law.
if sum((p.Results.eyePoseUB(1:3) - p.Results.eyePoseLB(1:3))==0) < 1
    warning('eyePoseEllipseFit:underconstrainedSearch','The eye pose search across possible eye rotations is underconstrained');
end

% If we have a glintCoord, make sure it is the right vector orientation
glintCoord = p.Results.glintCoord;
if ~isempty(glintCoord)
    if size(glintCoord,2)~=2
        error('eyePoseEllipseFit:glintCoordFormat','The optional glintCoord must be an nx2 vector');
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
unconstrainedEllipse = pupilEllipseFit([Xp,Yp]);

% Various degenerate sets of
% perimeter points can cause the ellipse fit to fail and return nans
if any(isnan(unconstrainedEllipse))
    % The fit failed. Sythesize an ellipse vector that has as its center
    % the mean of the X and Y positions of the perimeter points.
    unconstrainedEllipse = nan(1,5);
    unconstrainedEllipse(1:2) = [mean(Xp) mean(Yp)];
end


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
    % CoP.
    displacePix = (unconstrainedEllipse(1:2)-CoP);
    
    % If the ratio of the displacePix is extreme, or one or more of the
    % displacePix values is small, then the eye could be close to the
    % center of projection for one of the rotations. In this case, give the
    % same, unitary value to the displaceScaled variable
    if (displacePix(1) / displacePix(2)) > 2 || (displacePix(1) / displacePix(2)) < 0.5 || min(displacePix) < 1
        displaceScaled = sign(displacePix);
    else
        displaceScaled = sign(displacePix) .* (abs(displacePix) ./ max(abs(displacePix)));
    end
    
    % Handle the direction of eye rotation
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

% define some search options for fminbnd
options = optimset('fminsearch');
options.Display = 'off'; % Silencio
options.FunValCheck = 'off'; % The objective will return nans when the eyepose is not valid
options.TolX = p.Results.eyePoseTol; % eyePose search precision
options.MaxFunEvals = p.Results.eyePoseEllipseFitFunEvals; % This has the largest effect on search time

% Turn off warnings that can arise during the search
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% Clear the warning buffer
lastwarn('');

% Initialize the output variable
eyePose = x0;

% Search with a nested objective function
[eyePose(notLocked), RMSE] = fminsearch(@objFun, x0(notLocked), options);
    function fVal = objFun(x)
        % Combine the locked params with the searched params
        candidateEyePose = x0;
        candidateEyePose(notLocked) = x;
        % This is the value returned if no valid eye pose is found
        objFailValue = nan;
        % If any of the params exceed a boundary, return with the failed
        % objValue
        if any(candidateEyePose>eyePoseUB) || any(candidateEyePose<eyePoseLB)
            fVal = objFailValue;
            return
        end
        % Obtain the entrance pupil ellipse for this eyePose
        [candidateEllipse, candidateGlint] = projectModelEye(candidateEyePose, sceneGeometry);
        % Check for the case in which the transparentEllipse contains nan
        % values, which can arise if there were an insufficient number of
        % pupil border points remaining after refraction to define an
        % ellipse.
        if any(isnan(candidateEllipse))
            % Set fVal to something arbitrarily large
            fVal = objFailValue;
        else
            % This is the RMSE of the distance values of the boundary
            % points to the ellipse fit.
            explicitEllipse = ellipse_transparent2ex(candidateEllipse);
            if isempty(explicitEllipse)
                fVal = objFailValue;
            else
                if any(isnan(explicitEllipse))
                    fVal = objFailValue;
                else
                    fVal = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
                end
            end
        end
        % Check the match to the glint. The fit error is inflated by the
        % normed mismatch in pixels between the modeled and observed
        % glint(s).
        if ~isempty(glintCoord)
            if isempty(candidateGlint)
                fVal = objFailValue;
            else
                fVal = fVal * (1+norm(glintCoord - candidateGlint));
            end
        end
    end % local objective function

% Update the fittedEllipse with the solution parameters
fittedEllipse = projectModelEye(eyePose, sceneGeometry);

% Check if the fit is at a bound for any parameter that is not locked
fitAtBound = any([any(abs(eyePose(notLocked)-eyePoseLB(notLocked)) < p.Results.eyePoseTol) any(abs(eyePose(notLocked)-eyePoseUB(notLocked)) < p.Results.eyePoseTol)]);

% Restore the warning state
warning(warningState);

end % eyePoseEllipseFit


