function virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry, functionDirPath )
% Function handles to ray tracing equations
%
% Syntax:
%  virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry )
%
% Description:
%   This routine returns a handle to a function that is used to calculate
%   the location of a virtual image point that has passed through an
%   optical system. If the key-value 'functionDirPath' is set, then the
%   handle will be to a compiled mex function that has been written to
%   disk, otherwise, the handle will be to a matlab function in memory.
%   The former executes ~100x faster.
%
%   The virtualImageFunc is assembled step-wise from more elementary
%   algorithms:
%       traceOpticalSystem - 2D ray tracing through the cornea and any
%           corrective lenses
%       calcCameraNodeDistanceError2D - 2D distance of ray intersection on 
%           camera plane from camera node
%       calcVirtualImageRay - Returns the unit vector virtual image ray for
%           the initial depth position
%   This main routine calls out to local functions that assemble these
%   elementary components. Each elementary component is either saved as a
%   file at the specified location, or maintained as a function in memory.
%   Once these components are assembled, a handle is made to the function
%   virtualImageFunc, which uses these elementary components.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure. Critically, this
%                           includes an optical system.
%
% Optional key-value pairs:
%  'functionDirPath'     - Character string, default empty. If set this
%                           defines the location in which the compiled
%                           function is writen.
%  'cleanUpCompileDir'    - Logical, default true. If file path is
%                           provided, this flag determines if the
%                           intermediate compilation products are deleted.
%
% Outputs:
%   virtualImageFuncPointer - Structure. Includes the fields:
%                           'handle' - handle for the function.
%                           'path' -  full path to the stored mex file; set
%                               to empty if stored only in memory.
%                           'opticalSystem' - the optical system used to
%                               generate the function.
%
% Examples:
%{
    % Basic example with file caching of the functions
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, '/tmp/demo_virtualImageFunc' );
%}
%{
    % Demonstrate how the time it takes to perform the symbolic variable
    % calculations grows geometrically with the number of surfaces in the
    % optical system.

    % Obtain a default sceneGeometry. 
    sceneGeometry = createSceneGeometry();
    % Define the virtual image function
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(1)=toc;
    n(1)=size(sceneGeometry.opticalSystem,1);
    % Add a contact lens (one additional surface)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
    % Define the ray tracing functions 
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(2)=toc;
    n(2)=size(sceneGeometry.opticalSystem,1);
    % Add a spectacle lens (two additional surfaces)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',-2);
    % Define the ray tracing functions 
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(3)=toc;
    n(3)=size(sceneGeometry.opticalSystem,1);
    % Plot the timing results
    plot(n,t,'*r');
    xlabel('# of surfaces in optical model');
    ylabel('time to assemble ray tracing funcs [secs]');
%}




% Create a directory for the compiled functions
if ~isempty(functionDirPath)
    compileDir = functionDirPath;
    if ~exist(compileDir,'dir')
        mkdir(compileDir)
    end
else
    compileDir = [];
end


    syms z
    assume(z,'real');
    syms h
    assume(h,'real');
    syms theta
    assume(theta,'real');
    outputRayFunc = rayTraceCenteredSurfaces([z h], theta, sceneGeometry.opticalSystem);
    unity = 1; zero = 0;
    outputRayFunc = subs(outputRayFunc);
    functionFileName = fullfile(compileDir,'traceOpticalSystem');
    traceOpticalSystem = matlabFunction(outputRayFunc,'File',functionFileName,'Vars',{z h theta});    
    % Add saved function files to path
    addpath(compileDir,'-end');


% A functionDirPath has been defined, so we will compile the function as a
% mex file and place it on the path
    % Define some argument variables so that the compiler can deduce
    % variable types
    args = {[0, 0, 0], [0,0,0,0], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters};
    % Change to the compile directory
    initialDir = cd(compileDir);
    % Compile the mex file
    codegen -o virtualImageFuncMex virtualImageFuncPreMex -args args
    % Identify the compiled mex file, the suffix of which will vary
    % depending upon the operating system
    fileLocation = dir('virtualImageFuncMex.*');
    % Clean up the compile dir, if requested
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


end % compileVirtualImageFunc -- MAIN






