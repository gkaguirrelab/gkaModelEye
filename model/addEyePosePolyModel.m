function sceneGeometry = addEyePosePolyModel(sceneGeometry, varargin)
% Generate a polynomial model relating ellipse params to eyePose
%
% Syntax:
%  sceneGeometry = addEyePosePolyModel(sceneGeometry, varargin)
%
% Description:
%	This routine computes the pupil ellipse that is produce by a range of
%	eyePose values. The relationship between pupil ellipse values and
%	eyePose values are fit with a high-order polynomial model. This model
%	is then used to provide an initial x0 value for searches performed in
%	eyePoseEllipse fit and pupilProjection_inv.
%
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePoses for the search [azimuth,
%                           elevation, torsion, stop radius].
%  'gridDensity'          - Scalar. How densely the eyePose grid of values
%                           is sampled. Computation time will grow
%                           geometrically with grid density.
%  'estimateCompTime'     - Logical. If set to true, the routine will exit
%                           without performing the computation and will
%                           instead report the expected time to compute the
%                           grid for these parameter settings.
%
% Outputs:
%   sceneGeometry         - Structure. The input structure is copied over
%                           to the output, with the addition of the field
%                           polyModel.
%
% Examples:
%{
    % ETTBSkip -- This takes about 5 minutes to run
    % Generate the sceneGeometry and the polymodel
    sceneGeometry = createSceneGeometry();
    sceneGeometry = addEyePosePolyModel(sceneGeometry,'verbose',true);
    % Assume an eyePose and create a forward projection ellipse
    eyePose = [-10 5 0 2];
    pupilEllipse = pupilProjection_fwd(eyePose,sceneGeometry);
    % Given the pupilEllipse, reconstruct the eyePose using the polymodel
    x = zeros(1,4);
    x(1) = polyvaln(sceneGeometry.polyModel.azimuth,pupilEllipse);
    x(2) = polyvaln(sceneGeometry.polyModel.elevation,pupilEllipse);
    x(4) = polyvaln(sceneGeometry.polyModel.stopRadius,pupilEllipse);
    % Report the difference between the veridical and estimated eyePose
    eyePose - x

%}
%{
    % Estimate computation time
    sceneGeometry = createSceneGeometry();
    addEyePosePolyModel(sceneGeometry,'gridDensity',85,'estimateCompTime',true);
%}
%{
    % ETTBSkip -- This takes about 5 minutes to run
    %% Demonstrate acceleration of the inverse projection
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    sceneGeometry = addEyePosePolyModel(sceneGeometry,'verbose',true);
    % Generate ellipses for some randomly selected eye poses
    nPoses = 20;
    eyePoses=[(rand(nPoses,1)-0.5)*40, (rand(nPoses,1)-0.5)*20, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    for pp = 1:nPoses
    	ellipseParams(pp,:) = pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    fprintf('\nTime to compute inverse projection model with pre-computed polymodel (average over %d projections):\n',nPoses);
    recoveredEyePoses = []; RMSEvals = [];
    tic
    for pp = 1:nPoses
    	[recoveredEyePoses(pp,:),RMSEvals(pp), ~, ~, nSearches(pp)] = pupilProjection_inv(ellipseParams(pp,:),sceneGeometry,'repeatSearchThresh',0.5);
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
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('gridDensity',50,@(x)(isscalar(x) && x>=15));
p.addParameter('polyModelOrder',6,@isscalar);
p.addParameter('estimateCompTime',false,@islogical);
p.addParameter('verbose', false, @islogical);

% Parse and check the parameters
p.parse(sceneGeometry, varargin{:});

% Anonymous function to produce exponentially spaced values
expVal = 20;
expspace = @(mn,mx,n) (mx-mn)/expVal * (10.^(linspace(0, log10(expVal+1), n)) - 1) + mn;

% Generate the grid of values
aziVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(1)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(1),round(p.Results.gridDensity/2))]);
eleVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(2)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(2),round(p.Results.gridDensity/2))]);
stopVals = linspace(p.Results.eyePoseLB(4),p.Results.eyePoseUB(4),round(p.Results.gridDensity/4));

nComps = length(aziVals)*length(eleVals)*length(stopVals);

% Check if we are just to estimate the computation time
if p.Results.estimateCompTime || p.Results.verbose
    tic
    for ii=1:10
        pupilProjection_fwd([0 0 0 1], sceneGeometry, 'nStopPerimPoints', 5);
    end
    compTimeMins = toc*nComps/ii/60;
    fprintf('Estimated computation time for these parameters is %2.1f mins\n',compTimeMins);
end

if p.Results.estimateCompTime
    fprintf('Returning without performing the computation\n');
    return
end

% Loop through the eyePoses and calculate the forward projection
eyePoses = [];
pupilEllipses = [];
tic
idx=1;
for xx=1:length(aziVals)
    for yy=1:length(eleVals)
        for zz=1:length(stopVals)
            eyePose = [aziVals(xx) eleVals(yy) 0 stopVals(zz)];
            pupilEllipse = pupilProjection_fwd(eyePose, sceneGeometry, 'nStopPerimPoints', 5);
            if ~any(isnan(pupilEllipse))
                eyePoses(idx,:) = eyePose;
                pupilEllipses(idx,:) = pupilEllipse;
                idx=idx+1;
            end
        end
    end
end
compTimeMins = toc/60;

% report completion of prediction grid
if p.Results.verbose
    fprintf('Finished grid construction (%2.1f mins). Now fitting polynomial model\n',compTimeMins);
end

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:singularMatrix');

% Fit the polynomial model that relates pupilEllipse parameters to each of
% the eye pose parameters
sceneGeometry.polyModel.azimuth = ...
    polyfitn( pupilEllipses, eyePoses(:,1),p.Results.polyModelOrder);
sceneGeometry.polyModel.elevation = ...
    polyfitn( pupilEllipses, eyePoses(:,2),p.Results.polyModelOrder);
sceneGeometry.polyModel.stopRadius = ...
    polyfitn( pupilEllipses, eyePoses(:,4),p.Results.polyModelOrder);

% Restore the warning state
warning(warningState);

% Create an adjustment field and set the values to zero
sceneGeometry.polyModel.adjust = zeros(1,5);

% Store the meta data
sceneGeometry.polyModel.meta = p.Results;
sceneGeometry.polyModel.meta.compTimeMins = compTimeMins;

% report completion
if p.Results.verbose
    fprintf('Done\n\n');
end

end % function



