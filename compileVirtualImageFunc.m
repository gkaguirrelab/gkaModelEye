function virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry, functionDirPath )
% Function handles to ray tracing equations
%
% Syntax:
%  virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry )
%
% Description:
%   This routine produces a compiled mex file for virtualImageFunc, places
%   the function on the matlab path, and returns a handle to function.
%   The compiled version executes ~100x faster.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure. Critically, this
%                           includes an optical system.
%   functionDirPath       - Character vector. Specifies the location in 
%                           which the compiled function is writen.
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
    % Basic example with file caching of the functions
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
%}


%% Create a directory for the compiled functions
if ~isempty(functionDirPath)
    compileDir = functionDirPath;
    if ~exist(compileDir,'dir')
        mkdir(compileDir)
    end
else
    compileDir = [];
end

%% Save the traceOpticalSystem function
% Create a stand-alone 
syms z
%assume(z,'real');
syms h
%assume(h,'real');
syms theta
%assume(theta,'real');
outputRayFunc = rayTraceCenteredSurfaces([z h], theta, sceneGeometry.opticalSystem);
unity = 1; zero = 0;
outputRayFunc = subs(outputRayFunc);
functionFileName = fullfile(compileDir,'traceOpticalSystem');
traceOpticalSystem = matlabFunction(outputRayFunc,'File',functionFileName,'Vars',{z h theta});
% Add saved function files to path
addpath(compileDir,'-end');

%% Compile virtualImageFunc
% Define some argument variables so that the compiler can deduce
% variable types
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






