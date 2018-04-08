function virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry, varargin )
% Function handles to ray tracing equations
%
% Syntax:
%  virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry )
%
% Description:
%   This routine produces a compiled mex file for virtualImageFunc, saves
%   the file at a default (/tmp/demo_virtualImageFunc), and places the
%   function on the MATLAB path. If a second input argument is passed, this
%   is taken as the full path to where the compiled function should be
%   placed. The routine returns a structure that contains the handle to the
%   function.
%
%   Calls to the compiled virtualImageFuncMex execute roughly ~50x faster
%   than the native virtualImageFunc routine.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure. Critically, this
%                           includes an optical system.
%
% Optional inputs:
%   functionDirPath       - Character vector. Specifies the location in 
%                           which the compiled function is writen. Defaults
%                           to '/tmp/demo_virtualImageFunc'.
%
% Outputs:
%   virtualImageFuncPointer - Structure. Includes the fields:
%                           'handle' - handle for the function
%                           'path' -  full path to the stored mex file
%                           'opticalSystem' - the optical system used to
%                               generate the function
%
% Examples:
%{
    % Basic example with a compiled virtualImageFunc
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
    [virtualEyeWorldPoint, nodalPointIntersectError] = sceneGeometry.virtualImageFunc.handle( [sceneGeometry.eye.pupilCenter(1) 2 0], [0 0 0 2], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters )
%}
%{
    % Compare computation time for MATLAB and compiled C code
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
    nComputes = 1000;
    fprintf('\nTime to execute the virtualImageFunc (average over %d projections):\n',nComputes);
    % Native function
    tic
    for ii=1:nComputes
        virtualImageFunc( [sceneGeometry.eye.pupilCenter(1) 2 0], [0 0 0 2], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters );
    end
    msecPerCompute = toc / nComputes * 1000;
    fprintf('\tUsing the MATLAB function: %4.2f msecs.\n',msecPerCompute);
    % Compiled function
    tic
    for ii=1:nComputes
        sceneGeometry.virtualImageFunc.handle( [sceneGeometry.eye.pupilCenter(1) 2 0], [0 0 0 2], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters );
    end
    msecPerCompute = toc / nComputes * 1000;
    fprintf('\tUsing the compiled function: %4.2f msecs.\n',msecPerCompute);
%}

%% input parser
p = inputParser;

% Required
p.addRequired('sceneGeometry',@isstruct);
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
args = {[0,0,0], [0,0,0,0], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters};
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
addpath(compileDir,'-end');
% Change back to the initial directory
cd(initialDir);
% Return the path to the function as the output
virtualImageFuncPointer.handle = @virtualImageFuncMex;
virtualImageFuncPointer.path = fullfile(fileLocation.folder,fileLocation.name);
virtualImageFuncPointer.opticalSystem = sceneGeometry.opticalSystem;
% Save a copy of this variable in the function directory. The saved
% variable may be used to re-instantiate the function at a later point.
filePath = fullfile(fileLocation.folder,'virtualImageFuncPointer');
save(filePath,'virtualImageFuncPointer');

end % compileVirtualImageFunc






