function compileVirtualImageFunc( varargin )
% Compiles the virtualImageFunc and saves it to disk
%
% Syntax:
%  virtualImageFuncStruct = compileVirtualImageFunc( sceneGeometry, varargin )
%
% Description:
%   This routine produces a compiled mex file for virtualImageFunc, saves
%   the file at the specified disk location, and places the function on the
%   MATLAB path.
%
%   The default save location is within userpath(), at:
%
%       /toolboxes/transparentTrack/code/bin' 
%
%   Calls to the compiled virtualImageFuncMex execute roughly ~13x faster
%   than the native virtualImageFunc routine.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure.
%
% Optional key/value pairs:
%  'functionDirPath'      - Character vector. Specifies the location in 
%                           which the compiled function is writen.
%  'replaceExistingFunc'  - Logical, default false. If set to true, any
%                           existing versions of virtualImageFuncMex will
%                           be removed from the path and a new version will
%                           be created.
%
% Outputs:
%   none
%
% Examples:
%{
    % Basic example with a compiled virtualImageFunc
    compileVirtualImageFunc
    sceneGeometry = createSceneGeometry();
    % Assemble the args for the virtualImageFunc
    args = {sceneGeometry.cameraExtrinsic.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.opticalSystem.p1p2, ...
    	sceneGeometry.refraction.opticalSystem.p1p3};
    [virtualEyeWorldPoint, nodalPointIntersectError] = sceneGeometry.refraction.handle( [-3.7 2 0], [0 0 0 2], args{:} );
    % Test output against cached value
    virtualEyeWorldPointCached = [ -3.700000000000000   2.260297326476055  -0.000000817841759];
    assert(max(abs(virtualEyeWorldPoint - virtualEyeWorldPointCached)) < 1e-6)
%}
%{
    % Compare computation time for MATLAB and compiled C code
    % Turn off warning that we are replacing the compiled func
    nComputes = 100;
    fprintf('\nTime to execute virtualImageFunc (average over %d projections):\n',nComputes);
    % Native function
    sceneGeometry = createSceneGeometry('forceMATLABVirtualImageFunc',true);
    % Assemble the args for the virtualImageFunc
    args = {sceneGeometry.cameraExtrinsic.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.opticalSystem.p1p2, ...
    	sceneGeometry.refraction.opticalSystem.p1p3};
    tic
    for ii=1:nComputes
        virtualImageFunc( [-3.7 2 0], [0 0 0 2], args{:} );
    end
    msecPerComputeNative = toc / nComputes * 1000;
    fprintf('\tUsing the MATLAB function: %4.2f msecs.\n',msecPerComputeNative);
    % Compiled function
    sceneGeometry = createSceneGeometry();
    tic
    for ii=1:nComputes
        sceneGeometry.refraction.handle( [-3.7 2 0], [0 0 0 2], args{:} );
    end
    msecPerComputeCompile = toc / nComputes * 1000;
    fprintf('\tUsing the compiled function: %4.2f msecs.\n',msecPerComputeCompile);
%}

%% input parser
p = inputParser;

% Optional
p.addParameter('functionDirPath',fullfile(userpath(),'toolboxes','transparentTrack/code/bin'),@(x) ischar(x));
p.addParameter('replaceExistingFunc',false,@islogical);

% parse
p.parse(varargin{:});
functionDirPath = p.Results.functionDirPath;


%% Test if the function exists
% If we have not been asked to replace an existing function, test if the
% compiled function exists. If so, exit.
if ~p.Results.replaceExistingFunc
    if exist('virtualImageFuncMex')==3
        return
    end
end


%% Create a directory for the compiled functions
if ~isempty(functionDirPath)
    compileDir = functionDirPath;
    if ~exist(compileDir,'dir')
        mkdir(compileDir)
    end
else
    compileDir = [];
end


%% Remove pre-existing functions from the path
% Detect the case in which the current directory itself contains a compiled
% virtualImageFuncMex file, in which case the user needs to change
% directories
if strcmp(pwd(),fileparts(which('virtualImageFuncMex')))
    error('compileVirtualImageFunc:dirConflict','The current folder itself contains a compiled virtualImageFunc. Change directories to avoid function shadowing.')
end

% Remove any existing versions of the virtualImageFuncMex from the path.
notDoneFlag = true;
removalsCounter = 0;
tooManyRemovals = 4;
while notDoneFlag
    funcPath = which('virtualImageFuncMex');
    if isempty(funcPath)
        notDoneFlag = false;
    else
        warning('compileVirtualImageFunc:previousFunc','Removing a previous virtualImageFuncMex from the path');
        rmpath(fileparts(funcPath));
        removalsCounter = removalsCounter+1;
    end
    if removalsCounter == tooManyRemovals
        error('compileVirtualImageFunc:tooManyRemovals','Potentially stuck in a loop trying to remove previous virtualImageFuncMex functions from the path.')
    end
end


%% Compile virtualImageFunc
% Define argument variables so the compiler can deduce variable types
sceneGeometry = createSceneGeometry();
dynamicArgs = {[0,0,0], [0,0,0,0]};
staticArgs = {sceneGeometry.cameraExtrinsic.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.opticalSystem.p1p2, ...
    	sceneGeometry.refraction.opticalSystem.p1p3};

args = [dynamicArgs, staticArgs{:}];
% Change to the compile directory
initialDir = cd(compileDir);
% Compile the mex file
codegen -o virtualImageFuncMex virtualImageFunc -args args
% Identify the compiled mex file, the suffix of which will vary
% depending upon the operating system
fileLocation = dir('virtualImageFuncMex.*');
% Clean up the compile dir
rmdir('codegen', 's');
% Refresh the path to add the compiled function
addpath(compileDir,'-begin');
% Change back to the initial directory
cd(initialDir);


end % compileVirtualImageFunc






