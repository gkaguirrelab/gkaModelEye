function virtualImageFuncStruct = compileVirtualImageFunc( sceneGeometry, varargin )
% Compiles the virtualImageFunc and saves it to disk
%
% Syntax:
%  virtualImageFuncStruct = compileVirtualImageFunc( sceneGeometry, varargin )
%
% Description:
%   This routine produces a compiled mex file for virtualImageFunc, saves
%   the file at the specified disk location, and places the function on the
%   MATLAB path. If a second input argument is not passed, a default save
%   location is used (/tmp/demo_virtualImageFunc). The routine returns a
%   structure that contains the handle to the function.
%
%   Calls to the compiled virtualImageFuncMex execute roughly ~50x faster
%   than the native virtualImageFunc routine.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure.
%
% Optional inputs:
%   functionDirPath       - Character vector. Specifies the location in 
%                           which the compiled function is writen. Defaults
%                           to '/tmp/demo_virtualImageFunc'.
%
% Outputs:
%   virtualImageFuncStruct - Structure. Includes the fields:
%                           'handle' - handle for the function
%                           'path' -  full path to the stored mex file
%                           'opticalSystem' - the optical system(s) used to
%                               generate the function
%
% Examples:
%{
    % Basic example with a compiled virtualImageFunc
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
    [virtualEyeWorldPoint, nodalPointIntersectError] = sceneGeometry.virtualImageFunc.handle( [-3.7 2 0], [0 0 0 2], sceneGeometry.virtualImageFunc.args{:} );
    % Test output against value computed on April 10, 2018
    virtualEyeWorldPointStored = [-3.7000    2.2553    0.0000];
    assert(max(abs(virtualEyeWorldPoint - virtualEyeWorldPointStored)) < 1e-4)
%}
%{
    % Compare computation time for MATLAB and compiled C code
    sceneGeometry = createSceneGeometry();
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
    compileTime = toc;
    nComputes = 1000;
    clc
    fprintf('\nTime to execute virtualImageFunc (average over %d projections):\n',nComputes);
    % Native function
    tic
    for ii=1:nComputes
        virtualImageFunc( [-3.7 2 0], [0 0 0 2], sceneGeometry.virtualImageFunc.args{:} );
    end
    msecPerComputeNative = toc / nComputes * 1000;
    fprintf('\tUsing the MATLAB function: %4.2f msecs.\n',msecPerComputeNative);
    % Compiled function
    tic
    for ii=1:nComputes
        sceneGeometry.virtualImageFunc.handle( [-3.7 2 0], [0 0 0 2], sceneGeometry.virtualImageFunc.args{:} );
    end
    msecPerComputeCompile = toc / nComputes * 1000;
    fprintf('\tUsing the compiled function: %4.2f msecs.\n',msecPerComputeCompile);
    % Calculate the break-even point for perfoming the compilation
    breakEvenComputes = round(compileTime/((msecPerComputeNative-msecPerComputeCompile)/1000));
    fprintf('\nGiven a compile time of %4.2f seconds, break-even occurs when performing >%d computes.\n',compileTime,breakEvenComputes);
%}

%% input parser
p = inputParser;

% Required
p.addRequired('sceneGeometry',@isstruct);

% Optional
p.addOptional('functionDirPath',fullfile(filesep,'tmp','demo_virtualImageFunc'),@(x) ischar(x));

% parse
p.parse(sceneGeometry, varargin{:});
functionDirPath = p.Results.functionDirPath;


%% Create a directory for the compiled functions
if ~isempty(functionDirPath)
    compileDir = functionDirPath;
    if ~exist(compileDir,'dir')
        mkdir(compileDir)
    end
else
    compileDir = [];
end


%% Compile virtualImageFunc
% Define argument variables so the compiler can deduce variable types
args = {[0,0,0], [0,0,0,0], sceneGeometry.virtualImageFunc.args{:}};
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

% Return the path to the function as the output
virtualImageFuncStruct.handle = @virtualImageFuncMex;
virtualImageFuncStruct.path = fullfile(fileLocation.folder,fileLocation.name);
virtualImageFuncStruct.opticalSystem = sceneGeometry.virtualImageFunc.opticalSystem;

% Remake the args
virtualImageFuncStruct.args = {...
    sceneGeometry.cameraExtrinsic.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.virtualImageFunc.opticalSystem.p1p2, ...
    sceneGeometry.virtualImageFunc.opticalSystem.p1p3};

% Save a copy of this variable in the function directory. The saved
% variable may be used to re-instantiate the function at a later point.
filePath = fullfile(fileLocation.folder,'virtualImageFuncStruct');
save(filePath,'virtualImageFuncStruct');

end % compileVirtualImageFunc






