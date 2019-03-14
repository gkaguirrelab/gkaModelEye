function sceneGeometry = calcEyePoseGrid(sceneGeometry, varargin)
% Measure the pupil properties for a range of eye poses
%
% Syntax:
%  eyePoseGrid = calcEyePoseGrid(sceneGeometry, varargin)
%
% Description:
%	This routine computes and stores the pupil ellipse that is produce by a
%	range of eyePose values. These results are then stored in the field
%	sceneGeometry.eyePoseGrid. The grid is used to provide an initial x0
%	value for searches performed in eyePoseEllipse fit and
%	pupilProjection_inv.
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
%                           eyePoseGrid.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    sceneGeometry = calcEyePoseGrid(sceneGeometry);
%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('gridDensity',25,@isscalar);
p.addParameter('estimateCompTime',false,@islogical);

% Parse and check the parameters
p.parse(sceneGeometry, varargin{:});

% Anonymous function to produce exponentially spaced values
expVal = 20;
expspace = @(mn,mx,n) (mx-mn)/expVal * (10.^(linspace(0, log10(expVal+1), n)) - 1) + mn;

% Generate the grid of values
aziVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(1)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(1),round(p.Results.gridDensity/2))]);
eleVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(2)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(2),round(p.Results.gridDensity/2))]);
stopVals = linspace(p.Results.eyePoseLB(4),p.Results.eyePoseUB(4),round(p.Results.gridDensity/4));

% Check if we are just to estimate the computation time
if p.Results.estimateCompTime
    tic
    for ii=1:5
        pupilProjection_fwd([0 0 0 1], sceneGeometry, 'nStopPerimPoints', 5);
    end
    compTimeMins = toc*length(aziVals)*length(eleVals)*length(stopVals)/5/60;
    fprintf('Grid computation time for these parameters is %2.2f mins\n',compTimeMins);
    return
end

% Loop through the eyePoses and calculate the forward projection
eyePoses = [];
pupilEllipses = [];
tic
for xx=1:length(aziVals)
    for yy=1:length(eleVals)
        for zz=1:length(stopVals)
            eyePose = [aziVals(xx) eleVals(yy) 0 stopVals(zz)];
            pupilEllipse = pupilProjection_fwd(eyePose, sceneGeometry, 'nStopPerimPoints', 5);
            eyePoses(end+1,:)=eyePose;
            pupilEllipses(end+1,:)=pupilEllipse;
        end
    end
end
compTimeMins = toc/60;

% Store the results
sceneGeometry.eyePoseGrid.pupilEllipses = pupilEllipses;
sceneGeometry.eyePoseGrid.eyePoses = eyePoses;
sceneGeometry.eyePoseGrid.maxEllipseVals = max(abs(pupilEllipses),[],1);
sceneGeometry.eyePoseGrid.meta = p.Results;
sceneGeometry.eyePoseGrid.compTimeMins = compTimeMins;

end % function



