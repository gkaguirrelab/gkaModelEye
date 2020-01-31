function compileInverseRayTrace( varargin )
% Compiles the findPupilRay and findGlintRay functions
%
% Syntax:
%  compileInverseRayTrace
%
% Description:
%   This routine produces a compiled mex files for findPupilRay and
%   findGlintRay, saves the files at the specified disk location, and
%   places the functions on the MATLAB path.
%
%   The default save location is the directory that contains this function.
%
%   Calls to the compiled functions execute roughly ~30x faster than the
%   native routines.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'functionDirPath'      - Character vector. Specifies the location in 
%                           which the compiled function is writen.
%  'replaceExistingFunc'  - Logical, default false. If set to true, any
%                           existing versions of findPupilRayMex will
%                           be removed from the path and a new version will
%                           be created.
%
% Outputs:
%   none
%
% Examples:
%{
    % Confirm that compiled and native findPupilRay yield same value
    sceneGeometry = createSceneGeometry();
    % Assemble the args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.stopToCamera.opticalSystem};
    inverseRayNative = findPupilRay( [sceneGeometry.eye.stop.center(1) 2 0], [-5 10 0 2], args{:} );
    inverseRayCompiled = findPupilRayMex( [sceneGeometry.eye.stop.center(1) 2 0], [-5 10 0 2], args{:} );
    % Test if the outputs agree
    assert(max(max(abs(inverseRayNative - inverseRayCompiled))) < 1e-6)
%}
%{
    % Compare computation time for MATLAB and compiled code
    nComputes = 100;
    fprintf('\nTime to execute findPupilRay (average over %d projections):\n',nComputes);
    % Native function
    sceneGeometry = createSceneGeometry();
    % Assemble the args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.stopToCamera.opticalSystem};
    % Native matlab function
    tic
    for ii=1:nComputes
        findPupilRay( [-3.7 2 0], [0 0 0 2], args{:} );
    end
    msecPerComputeNative = toc / nComputes * 1000;
    fprintf('\tUsing the MATLAB function: %4.2f msecs.\n',msecPerComputeNative);
    % Compiled MEX function
    tic
    for ii=1:nComputes
        findPupilRayMex( [-3.7 2 0], [0 0 0 2], args{:} );
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
    if exist('findPupilRayMex')==3
        warning('compileInverseRayTrace:functionExists','findPupilRayMex already exists; set replaceExistingFunc to true to over-write.')
        return
    end
end


%% Error if the function dir does not exist
if ~exist(functionDirPath,'dir')
    error('compileInverseRayTrace:dirDoesNotExist','The specified function directory does not exist.')
end


%% Remove pre-existing functions from the path
% Detect the case in which the current directory itself contains a compiled
% findPupilRayMex file, in which case the user needs to change
% directories
if strcmp(pwd(),fileparts(which('findPupilRayMex')))
    error('compileInverseRayTrace:dirConflict','The current folder itself contains a compiled findPupilRay. Change directories to avoid function shadowing.')
end

% Remove any existing versions of the findPupilRayMex from the path.
notDoneFlag = true;
removalsCounter = 0;
tooManyRemovals = 4;
while notDoneFlag
    funcPath = which('findPupilRayMex');
    if isempty(funcPath)
        notDoneFlag = false;
    else
        warning('compileInverseRayTrace:previousFunc','Removing a previous findPupilRayMex from the path');
        rmpath(fileparts(funcPath));
        removalsCounter = removalsCounter+1;
    end
    if removalsCounter == tooManyRemovals
        error('compileInverseRayTrace:tooManyRemovals','Potentially stuck in a loop trying to remove previous findPupilRayMex functions from the path.')
    end
end


%% Define argument variables
% This is so the compiler can deduce variable types

% Create a sceneGeometry.
sceneGeometry = createSceneGeometry();
% Define the form of the dynamicArgs (the eyePoint and the eyePose)
dynamicArgs = {[0,0,0], [0,0,0,0]};
% Define the form of the staticArgs (which are sceneGeometry components)
staticArgs = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.stopToCamera.opticalSystem};
% Assemble the full args
args = [dynamicArgs, staticArgs{:}];


%% Compile and clean up
% Change to the compile directory
initialDir = cd(functionDirPath);
% Compile the mex file
codegen -o findPupilRayMex findPupilRay -args args
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


end % compileInverseRayTrace






