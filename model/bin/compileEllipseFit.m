function compileEllipseFit( varargin )
% Compiles the ellipseFitRobust function
%
% Syntax:
%  compileEllipseFit
%
% Description:
%   This routine produces a compiled MEX function for fitting an ellipse to
%   a set of points. This is saved at the specified disk location, and
%   placed on the MATLAB path.
%
%   The default save location is the directory that contains this function.
%
%   Calls to the compiled functions execute ~30x faster than the native
%   routines.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'functionDirPath'      - Character vector. Specifies the location where
%                           the compiled functions are saved.
%  'replaceExistingFunc'  - Logical, default false. If set to true, any
%                           existing version of the compiled functions for
%                           the current operating system will be removed
%                           from the path and a new version will be
%                           created.
%
% Outputs:
%   none
%
% Examples:
%{
    % Confirm that compiled and native ellipseFit_robust yield the same value
    eyePose = [10 -5 0 2.5];
    sceneGeometry=createSceneGeometry();
    [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );
    ellipseFitNative = pupilEllipseFit([Xp,Yp],@ellipsefit_robust);
    ellipseFitCompiled = pupilEllipseFit([Xp,Yp],@ellipsefit_robustMex);
    % Test if the outputs agree
    assert(max(max(abs(ellipseFitNative - ellipseFitCompiled))) < 1e-6)
%}
%{
    % Compare computation time for MATLAB and compiled code
    nComputes = 100;
    fprintf('Time to execute ellipseFit_robust (average over %d executions):\n',nComputes);
    % Assemble the inputs
    eyePose = [10 -5 0 2.5];
    sceneGeometry=createSceneGeometry();
    [ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry);
    [ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );
    ellipseFitCompiled = pupilEllipseFit([Xp,Yp],@ellipsefit_robustMex);
    % Native function
    tic
    for ii=1:nComputes
        pupilEllipseFit([Xp,Yp],@ellipsefit_robust);
    end
    msecPerComputeNative = toc / nComputes * 1000;
    fprintf('\tUsing the MATLAB function: %4.2f msecs.\n',msecPerComputeNative);
    % Compiled MEX function
    tic
    for ii=1:nComputes
        pupilEllipseFit([Xp,Yp],@ellipsefit_robustMex);
    end
    msecPerComputeCompile = toc / nComputes * 1000;
    fprintf('\tUsing the compiled function: %4.2f msecs.\n',msecPerComputeCompile);
%}


%% input parser
p = inputParser;

% Optional
p.addParameter('functionDirPath',fileparts(mfilename('fullpath')),@(x) ischar(x));
p.addParameter('replaceExistingFunc',false,@islogical);

% parse
p.parse(varargin{:});
functionDirPath = p.Results.functionDirPath;


%% Test if the function exists
% If we have not been asked to replace an existing function, test if the
% compiled function exists. If so, exit.
if ~p.Results.replaceExistingFunc
    % Exist returns 3 for 'MEX-file on your MATLAB search path'
    if exist('ellipsefit_robustMex','file')==3
        warning('compileEllipseFit:functionExists','ellipsefit_robustMex already exists; set replaceExistingFunc to true to over-write.')
        return
    end
end


%% Error if the function dir does not exist
if ~exist(functionDirPath,'dir')
    error('compileEllipseFit:dirDoesNotExist','The specified function directory does not exist.')
end


%% Remove pre-existing functions from the path
% Detect the case in which the current directory itself contains a compiled
% findPupilRayMex file, in which case the user needs to switch directories
if strcmp(pwd(),fileparts(which('ellipsefit_robustMex')))
    error('compileEllipseFit:dirConflict','The current folder itself contains a compiled ellipsefit_robust. Change directories to avoid function shadowing.')
end

% Remove any existing versions of the compiled code from the path.
notDoneFlag = true;
removalsCounter = 0;
tooManyRemovals = 4;
while notDoneFlag
    funcPath = which('ellipsefit_robustMex');
    if isempty(funcPath)
        notDoneFlag = false;
    else
        warning('compileEllipseFit:previousFunc','Removing previous compiled functions from the path');
        rmpath(fileparts(funcPath));
        removalsCounter = removalsCounter+1;
    end
    if removalsCounter == tooManyRemovals
        error('compileEllipseFit:tooManyRemovals','Potentially stuck in a loop trying to remove previous compiled functions from the path.')
    end
end


%% Define argument variables
% This is so the compiler can deduce variable types
args = {zeros(6,6),zeros(6,6)};

%% Compile and clean up
% Change to the compile directory
initialDir = cd(functionDirPath);
% Compile the mex file
codegen -o ellipsefit_robustMex ellipsefit_robust -args args
% Clean up the compile dir. Turn off warnings regarding the removal of
% these files
warnState = warning();
warning('Off','MATLAB:RMDIR:RemovedFromPath');
rmdir('codegen', 's');
warning(warnState);
% Refresh the path to add the compiled function
addpath(functionDirPath,'-begin');
% Change back to the initial directory
cd(initialDir);


end % compileEllipseFit






