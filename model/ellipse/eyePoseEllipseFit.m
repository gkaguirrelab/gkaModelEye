function [eyePose, cameraTrans, RMSE, fittedEllipse, fitAtBound, searchOutput, xHist] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, options)
% Return the eye pose that best fits an ellipse to the pupil perimeter
%
% Syntax:
%  [eyePose, cameraTrans, RMSE, fittedEllipse, fitAtBound, searchOutput, xHist] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry)
%
% Description:
%   The routine finds the eyePose parameters (azimuth, elevation, torsion,
%   stop radius) that produces a pupil ellipse that best fits a set of
%   pupil perimeter points.
%
%   If one or more glint coordinates are supplied, the eyePose used to fit
%   the pupil perimeter is constrained to also produce a modeled location
%   of the glint(s) that matches the supplied coordinates.
%
%   If a glint is provided, a search across translation of the camera is
%   also supported. If there is only one glint, then only camera in-plane
%   translation can be estimated, although not camera depth. If there is
%   more than one glint, then 3D camera translation will be estimated. The
%   routine will provide a warning if the bounds on the cameraTranslation
%   search are not appropriate for the glint information provided.
%
%   The search is constrained by the upper and lower bounds of the eyePose.
%   The default values specified here represent the physical boundaries of
%   the rotation model. Tighter, biologically informed constraints may be
%   passed by the calling function.
%
%   The rotation values returned in the eyePose variable are relative to
%   the primary position of the eye.
%
% Inputs:
%   Xp, Yp                - mx1 vectors of points to be fit. These provide
%                           the X and Y screen coordinates of the m points
%                           on the boundary of the pupil.
%   glintCoord            - A nx2 vector with the image coordinates of the
%                           n glints.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'eyePoseX0'            - A 1x4 vector that provides starting points for
%                           the search for the eyePose. If not defined, the
%                           starting point will be estimated.
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePose [azimuth, elevation,
%                           torsion, stop radius]. The default values here
%                           represent the physical limits of the projection
%                           model for azimuth, elevation, and stop radius.
%                           Torsion is constrained to zero by default.
%  'cameraTransX0'        - A 3x1 vector that provides starting points for
%                           camera translation, in units of mm. If  not
%                           defined, the starting point will be estimated.
%  'cameraTransBounds'    - A 3x1 vector that provides the symmetric bounds
%                           around x0 on camera translation [horizontal;
%                           vertical; depth]. Note that this manner of
%                           specifying bounds differs from what is used for
%                           the eyePose.
%  'eyePoseEllipseFitFunEvals' - Scalar. The maximum number of evals for
%                           the fmincon operation. The example below
%                           examines the trade-off between execution time
%                           and function accuracy for this parameter. I
%                           find that 250 evals takes 300ms per eyePose
%                           calculation on a laptop, and produces excellent
%                           accuracy with respect to the precision in
%                           empirical data.
%  'eyePoseTol'           - Scalar. The eyePose values will be searched to
%                           within this level of precision. Doesn't have
%                           too much of an effect upon execution time.
%  'glintTol'             - Scalar. This sets the tolerance (in units of
%                           pixels) on the non-linear constraint that
%                           forces the modeled and observed glint to match.
%  'boundTol'             - 1x7 vector. If a eyePose/cameraTrans value is
%                           within this distance of the bound, it is
%                           considered to have hit the bound. This is
%                           required as fmincon seems to avoid the bounds
%                           by a bit.
%
% Outputs:
%   eyePose               - A 1x4 vector with values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           degrees, and stop radius is in mm.
%   cameraTrans           - A 3x1 vector with values for [horizontal;
%                           vertical; depth] in units of mm. The values are
%                           relative to:
%                               sceneGeometry.cameraPosition.translation
%   RMSE                  - Root mean squared error of the distance of
%                           boundary points in the image to the fitted
%                           ellipse
%   fittedEllipse         - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   fitAtBound            - Logical. Indicates if any of the returned
%                           parameters (eye or camera) are at the upper or
%                           lower boundary.
%   searchOutput          - Structure. The output of the fmincon search.
%   xHist                 - A nx7 array, where n is equal to the number of
%                           objective function evaluations in the fmincon
%                           search. The values correspond to the eyePose
%                           and cameraTranslation examined in each search
%                           iteration.
%
% Examples:
%{
    % Basic example of recovering a simulated eyePose under the assumption
    % of no camera translation
    eyePose = [10 -5 0 2.5];
    sceneGeometry=createSceneGeometry();
    [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );
    eyePoseRecovered = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, 'cameraTransBounds', [0; 0; 0]);
    assert(max(abs(eyePose - eyePoseRecovered)) < 1e-2);
%}
%{
    % Eye pose in the setting of relative camera motion
    sceneGeometry=createSceneGeometry();
    eyePose = [10 0.5 0 2.5];
    cameraTrans = [3; -1; 0];
    [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry,'cameraTrans',cameraTrans);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );
    [eyePoseRecovered, cameraTransRecovered, RMSE, ~, ~, searchOutput] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry);
    assert(max(abs(eyePose-eyePoseRecovered)) < 1e-1);
    assert(max(abs(cameraTrans-cameraTransRecovered)) < 1e-1);
%}
%{
    % Explore the trade-off between accuracy and run time in the setting
    % of a relative camera translation
    sceneGeometry=createSceneGeometry();
    cameraTrans = [3; -1; 0];
    executionTime = []; errors = [];
    poseCount = 0;
    evalNum = 100:50:300;
    for azi = [-10 -.5 0.5 10]; for ele = [-10 -.5 0.5 10]
        eyePose = [azi ele 0 2.5];
        poseCount = poseCount+1;
        [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry,'cameraTrans',cameraTrans);
        [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );
        for ee = 1:length(evalNum)
            tic
            recoveredEyePose = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry,...
                'eyePoseEllipseFitFunEvals', evalNum(ee));
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
    plot(evalNum,timeByNumEvals,'-r','lineWidth',2)
    ylabel('Execution time per pose [secs]');
    yyaxis left
    plot(evalNum,errorByNumEvals,'-b','lineWidth',2)
    ylabel('Eye pose error [deg]');
    xlabel('Number of search evaluations');
    title('Performance of eyePoseEllipseFit across max fun evals');
%}

arguments
    Xp double
    Yp double
    glintCoord double
    sceneGeometry struct
    options.eyePoseX0 double = []
    options.eyePoseLB double = [-89, -89, 0, 0.1]
    options.eyePoseUB double = [89, 89, 0, 4]
    options.cameraTransX0 double = []
    options.cameraTransBounds double = [5; 5; 0]
    options.eyePoseEllipseFitFunEvals double {mustBeScalarOrEmpty} = 250
    options.eyePoseTol double {mustBeScalarOrEmpty} = 1e-3
    options.glintTol double {mustBeScalarOrEmpty} = 1
    options.boundTol double {mustBeVector} = [0.1 0.1 0.1 0.05 0.1 0.1 0.1]
end


%% Prepare return variables and check inputs

% Initialize the return variables
eyePose = [nan nan nan nan];
RMSE = nan;
fittedEllipse = [nan nan nan nan nan];
fitAtBound = false;

% Issue a warning if the bounds do not fully constrain at least one eye
% rotation parameter. This is because there are multiple combinations of
% the three axis rotations that can bring an eye to a destination.
% Typically, the torsion will be constrained with upper and lower bounds of
% zero (SEE: project/stages/addPseudoTorsion.m).
if sum((options.eyePoseUB(1:3) - options.eyePoseLB(1:3))==0) < 1
    warning('eyePoseEllipseFit:underconstrainedSearch','The eye pose search across possible eye rotations is underconstrained');
end

% Convert a nan glintCoord to empty
if any(isnan(glintCoord))
    glintCoord = [];
end

% If we have a glintCoord, make sure it is the right vector orientation
if ~isempty(glintCoord)
    if size(glintCoord,2)~=2
        error('eyePoseEllipseFit:glintCoordFormat','The glintCoord must be empty or an nx2 vector');
    end
end

% Clear the case of nans in the input
if any(isnan(options.eyePoseLB)) || ...
        any(isnan(options.eyePoseUB)) || ...
        any(any(isnan(options.eyePoseX0))) || ...
        any(isnan(options.cameraTransBounds)) || ...
        any(any(isnan(options.cameraTransX0)))
    return
end

% Clear the case of an empty perimeter
if isempty(Xp) || isempty(Yp)
    return
end

% Issue a warning if there are non-zero bounds on camera translation, but
% there is no glint. Without a glint, we can't estimate camera translation.
if isempty(glintCoord) && any(abs(options.cameraTransBounds) > 0)
    warning('eyePoseEllipseFit:underconstrainedSearch','No glint provided; cameraTrans search is under-constrained');
end

% Issue a warning if there are non-zero bounds on camera depth translation,
% but we only have one glint. Without two or more glints, we don't have
% much traction on depth.
if size(glintCoord,1)<2 && abs(options.cameraTransBounds(3)) > 0
    warning('eyePoseEllipseFit:underconstrainedSearch','Only one glint provided; cameraTrans depth search is under-constrained');
end


%% Define cameraTransX0 and bounds
if ~isempty(options.cameraTransX0)
    cameraTransX0 = options.cameraTransX0;
    cameraTransLB = options.cameraTransX0 - options.cameraTransBounds;
    cameraTransUB = options.cameraTransX0 + options.cameraTransBounds;
else
    
    % The bounds will be set about zero
    cameraTransLB = -options.cameraTransBounds;
    cameraTransUB = options.cameraTransBounds;
    
    % If we were given an eyePoseX0, use this as the pose of the eye for
    % the calculation. Otherwise, assume the eye is in primary position.
    if ~isempty(options.eyePoseX0)
        refPose = options.eyePoseX0;
    else
        refPose = [0 0 0 2];
    end
    
    % Probe to find the translation in pixels produced
    [~, probeCoord0] = projectModelEye(refPose, sceneGeometry, 'cameraTrans', [0;0;0]);
    [~, probeCoord1] = projectModelEye(refPose, sceneGeometry, 'cameraTrans', [1;0;0]);
    
    % If the glint was undefined, set the cameraTransX0 to zero
    if isempty(probeCoord0) || isempty(probeCoord1)
        
        cameraTransX0 = [0;0;0];
        
    else
        
        pixPerMm = mean(probeCoord0(:,1)-probeCoord1(:,1));
        cameraTransX0 = [ ...
            -(glintCoord(1,1)-probeCoord0(1,1))/pixPerMm; ...
            (glintCoord(1,2)-probeCoord0(1,2))/pixPerMm; ...
            0];
        
        % If there are two or more glints, determine the depth translation
        if size(glintCoord,1)>1
            [~, probeCoord2] = projectModelEye(refPose, sceneGeometry, 'cameraTrans', [0;0;1]);
            a = diff(probeCoord0(1:2,1))/2;
            b = diff(probeCoord2(1:2,1))/2;
            t = (diff(glintCoord(1:2,1))/2 - a)/(b-a);
            cameraTransX0(3) = t;
        end
        
    end
    
    % Make sure that X0 is in bounds
    cameraTransX0 = max([cameraTransLB cameraTransX0],[],2);
    cameraTransX0 = min([cameraTransUB cameraTransX0],[],2);
    
end



%% Define eyePoseX0 and bounds

% bounds
eyePoseLB = options.eyePoseLB;
eyePoseUB = options.eyePoseUB;

% Check if eyePoseX0 was passed
if isempty(options.eyePoseX0)
    
    % Define eyePoseX0
    eyePoseX0 = zeros(1,4);
    
    % Obtain the unconstrained ellipse fit to the perimeter to identify the
    % center of the entrance pupil.
    unconstrainedEllipse = pupilEllipseFit([Xp,Yp]);
    
    % If the fit failed, sythesize an ellipse vector that has as its center
    % the mean of the X and Y positions of the perimeter points.
    if any(isnan(unconstrainedEllipse))
        unconstrainedEllipse = nan(1,5);
        unconstrainedEllipse(1:2) = [mean(Xp) mean(Yp)];
    end

    % Construct an eyePoseX0 guess by probing the forward model. First
    % identify the center of projection
    rotationCenterEllipse = projectModelEye([0 0 0 2], sceneGeometry, 'cameraTrans',cameraTransX0);
    CoP = [rotationCenterEllipse(1),rotationCenterEllipse(2)];
    
    % Now the number of pixels that the pupil center is displaced from the
    % center of projection.
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
    probeEllipse = projectModelEye([displaceScaled 0 2],sceneGeometry,'cameraTrans',cameraTransX0);
    pixelsPerDeg = abs(probeEllipse(1:2)-CoP)./abs(displaceScaled);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Need to make the
    % second value negative to match the direction of rotation convention
    % for elevation.
    eyePoseX0(1) = displacePix(1)/pixelsPerDeg(1);
    eyePoseX0(2) = -displacePix(2)/pixelsPerDeg(2);
    
    % Force the angles within bounds
    eyePoseX0=[min([eyePoseUB(1:3); eyePoseX0(1:3)]) 0];
    eyePoseX0=[max([eyePoseLB(1:3); eyePoseX0(1:3)]) 0];
    
    % Estimate the pupil radius in pixels
    pupilRadiusPixels = max([abs(max(Xp)-min(Xp)) abs(max(Yp)-min(Yp))])/2;
    
    % Probe the forward model at the estimated pose angles to estimate the
    % pupil radius.
    probeEllipse=projectModelEye([eyePoseX0(1) eyePoseX0(2) eyePoseX0(3) 2], sceneGeometry,'cameraTrans',cameraTransX0);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    eyePoseX0(4) = pupilRadiusPixels/pixelsPerMM;
    
    % Ensure that eyePoseX0 lies within the bounds with a bit of headroom
    % so that the solver does not get stuck up against a bound.
    boundHeadroom = (eyePoseUB - eyePoseLB)*0.001;
    eyePoseX0=min([eyePoseUB-boundHeadroom; eyePoseX0]);
    eyePoseX0=max([eyePoseLB+boundHeadroom; eyePoseX0]);
    
    % Force any locked parameter to have the locked value
    eyePoseX0(eyePoseLB == eyePoseUB) = eyePoseUB(eyePoseLB == eyePoseUB);
    
else
    eyePoseX0 = options.eyePoseX0;
end


%% Set up variables and functions for the search

% define some search options for fmincon
opt_fmincon = optimset('fmincon');
opt_fmincon.Display = 'off'; % Silencio
opt_fmincon.FunValCheck = 'off'; % Tolerate nans from the objective
opt_fmincon.TolX = options.eyePoseTol; % eyePose search precision
opt_fmincon.TolCon = options.glintTol; % glint match precision
opt_fmincon.MaxFunEvals = options.eyePoseEllipseFitFunEvals; % This has the largest effect on search time

% Turn off warnings that can arise during the search
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% Clear the warning buffer
lastwarn('');

% Assemble the combined eyePose and cameraTrans x0 and bounds
x0 = [eyePoseX0 cameraTransX0'];
lb = [eyePoseLB cameraTransLB'];
ub = [eyePoseUB cameraTransUB'];

% Initialize nested variables
xHist = x0;
fHist = nan;
ceqHist = nan;
xLast = [];
cLast = [];
ceqLast = [];
fValLast = [];
ceqFirstFlag = true;

% Initialize an anonymous function for the full objective
fullObj = @(x) fullObjective(x,Xp,Yp,glintCoord,sceneGeometry);

% Search with a nested objective function
[x, RMSE, exitFlag, searchOutput, lambda] = fmincon(@objFun, x0,[],[],[],[],lb,ub,@nonlcon,opt_fmincon);

    % The fVal is the RMSE of the pupil ellipse to the pupil perimeter
    function fVal = objFun(x)
        if ~isequal(x,xLast)
            [fValLast,cLast,ceqLast] = fullObj(x);
            xLast = x;
            xHist(end+1,:) = x;
            fHist(end+1) = fValLast;
            ceqHist(end+1) = ceqLast;
        end
        fVal = fValLast;
        % This is some business to prevent fmincon from considering as
        % acceptable a ceq that is low relative to the x0 ceq
        if ceqFirstFlag
            ceqFirstFlag = false;
        else
            if ceqLast == 1e3
                ceqLast = 1e6;
            end
        end
    end

    % This is the mismatch of the observed and modeled glint
    function [c, ceq] = nonlcon(x)
        if ~isequal(x,xLast)
            [fValLast,cLast,ceqLast] = fullObj(x);
            xLast = x;
            xHist(end+1,:) = x;
            fHist(end+1) = fValLast;
            ceqHist(end+1) = ceqLast;
        end
        c = cLast;
        % This is some business to prevent fmincon from considering as
        % acceptable a ceq that is low relative to the x0 ceq
        ceq = ceqLast;
        if ceqFirstFlag
            ceqFirstFlag = false;
        else
            if ceq == 1e3
                ceq = 1e6;
            end
        end
    end

% If the search hit the max fun evals, look in the search history to find
% the best solution. We use the Lagrange multiplier reported at the end of
% the search to weight the non-linear constraint.
if exitFlag == 0
    [~,idx] = min(fHist + ceqHist*lambda.eqnonlin);
    RMSE = fHist(idx);
    ceq = ceqHist(idx);
    x = xHist(idx,:);
    fHist(end+1) = RMSE;
    ceqHist(end+1) = ceq;
    xHist(end+1,:) = x;
end

% Unpack the x parameters into the eyePose and cameraTrans
eyePose = x(1:4);
cameraTrans = x(5:7)';

% Update the fittedEllipse with the solution parameters
fittedEllipse = projectModelEye(eyePose, sceneGeometry, 'cameraTrans', cameraTrans);

% Check if the fit is within boundTol of a bound for any non-locked
% parameter.
notLocked = lb ~= ub;
fitAtBound = any([any(abs(x(notLocked)-lb(notLocked)) < options.boundTol(notLocked)) any(abs(x(notLocked)-ub(notLocked)) < options.boundTol(notLocked))]);

% Restore the warning state
warning(warningState);

end % eyePoseEllipseFit


%% LOCAL FUNCTION
function [fVal,c,ceq] = fullObjective(x,Xp,Yp,glintCoord,sceneGeometry)
% This function returns both the objective and non-linear constraint 

% Set the default values
fVal = 1e6;
c = [];
ceq = 0;

% Obtain the pupil ellipse for this eyePose
[candidateEllipse, candidateGlint] = ...
    projectModelEye(x(1:4), sceneGeometry, 'cameraTrans', x(5:7)');

% Check for the case in which the transparentEllipse contains nan values,
% which can arise if there were an insufficient number of pupil border
% points remaining after refraction to define an ellipse
if any(isnan(candidateEllipse))
    return
end
    
% Calculate the RMSE of the distance values of the boundary points to the
% ellipse fit. First cast the ellipse in explicit form.
explicitEllipse = ellipse_transparent2ex(candidateEllipse);

% Detect failures in the ellipse fit, which will cause the explicit ellipse
% parameters to be empty or contain nans.
if isempty(explicitEllipse)
    return
end
if any(isnan(explicitEllipse))
    return
end

% Obtain the fVal, which is the RMSE of the ellipse fit to the pupil
% perimeter points.
fVal = sqrt(mean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2,'omitmissing'));

% Now compute the constraint. If there is no glint, then we do not use the
% constraint.
if isempty(glintCoord)
    return
end

% Detect cases in which we were unable to obtain a valid glint
if isempty(candidateGlint)
    ceq = 1e3;
    return
end

% The non-linear constraint is the number of pixels by which model misses
% the glint location(s)
ceq = norm(glintCoord - candidateGlint);

end

