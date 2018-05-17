function [eyePose, RMSE] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Syntax:
%  [eyePose, RMSE] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs)
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
%
% Outputs:
%   eyePose               - A 1x4 matrix containing the best fitting eye
%                           parameters (azimuth, elevation, torsion, pupil
%                           radius)
%   RMSE                  - Root mean squared error of the distance of
%                           boundary point in the image to the fitted
%                           ellipse
%


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

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


%% Set bounds and x0
% Identify the center of projection.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

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
    % estimate the pupil radius. Here we do need ray tracing as it
    % has a substantial influence upon the area of the ellipse.
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


% Define an anonymous function for the objective
myObj = @(x) objfun(x, Xp, Yp, sceneGeometry);

% define some search options
options = optimoptions(@fmincon,...
    'Display','off');

% Perform the non-linear search .We sometimes obtain a singular matrix
% warning here; turn it off temporarily
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

[eyePose, RMSE, exitFlag] = ...
    fmincon(myObj, x0, [], [], [], [], eyePoseLB, eyePoseUB, [], options);

% If exitFlag==2, we might be in a local minimum; try again starting from
% a position close to the point found by the prior search
if exitFlag == 2
    x0 = eyePose+[1e-6 1e-6 0 1e-6];
    [eyePose, RMSE] = ...
        fmincon(myObj, x0, [], [], [], [], p.Results.eyePoseLB, p.Results.eyePoseUB, [], options);
end

% Restore the warning state
warning(warningState);


end % eyeParamEllipseFit


%% LOCAL FUNCTIONS
function fVal = objfun(x, Xp,Yp, sceneGeometry)
% Define the objective function
transparentEllipse = pupilProjection_fwd(x, sceneGeometry);
% Check for the case in which the transparentEllipse contains NAN values,
% which can happen when the eye pose is such that the border of the pupil
% would not be visible through the cornea. In this case, we return a
% realMax value for the fVal.
if any(isnan(transparentEllipse))
    fVal = realmax;
else
    % This is the RMSE of the distance values of the boundary points to
    % the ellipse fit.
    explicitEllipse = ellipse_transparent2ex(transparentEllipse);
    fVal = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
end
end % local objective function


