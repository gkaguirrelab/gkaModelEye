function sceneGeometry = loadSceneGeometry(sceneGeometryFileName, verbose)
% Load sceneGeometry and instantiate the ray tracing function
%
% Syntax:
%  sceneGeometry = loadSceneGeometry(sceneGeometryFileName, verbosity)
%
% Description:
%   Loads the sceneGeometry file from file and ensures that a valid
%   virtual image ray-tracing function is available.
%
% Inputs:
%   sceneGeometryFileName - Full path to the sceneGeometry file. If left
%                           empty, then an empty variable will be returned.
%   verbose               - Boolean. Defaults false if not passed.
%
% Outputs:
%   sceneGeometry         - The sceneGeometry structure
%

% Set the verbose flag to false if not passed
if nargin==1
    verbose = false;
end


if isempty(sceneGeometryFileName)
    sceneGeometry=[];
else
    % load the sceneGeometry structure
    dataLoad=load(sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    % instantiate the ray-tracing function
    if ~isempty(sceneGeometry.refraction)
        % The field is not empty, so we should have a ray tracing function.
        try
            funResult = functions(sceneGeometry.refraction.handle);
        catch
            error('Invalid function definition in sceneGeometry.refraction');
        end
        switch funResult.type
            case 'anonymous'
                % No further action needed
                if verbose
                    fprintf('Anonymous ray trace function available\n');
                end
            case 'simple'
                % Determine if the function exists
                if exist(func2str(sceneGeometry.refraction.handle))==0
                    % We need to add the function back to the path
                    addpath(fileparts(sceneGeometry.refraction.path),'-BEGIN')
                    % Check to make sure that it is now available
                    if exist(func2str(sceneGeometry.refraction.handle))==0
                        error('Unable to re-instantiate the ray tracing function')
                    end
                else
                    % If we have a compiled function, make sure that it is
                    % the right compiled function.
                    if exist(func2str(sceneGeometry.refraction.handle))==3
                        if ~strcmp(which(func2str(sceneGeometry.refraction.handle)), sceneGeometry.refraction.path)
                            % Attempt to remove the currently prioritized
                            % function. Silence a warning about it not
                            % existing, which can occur.
                            warnState=warning();
                            warning('off','MATLAB:rmpath:DirNotFound');
                            rmpath(fileparts(which(func2str(sceneGeometry.refraction.handle))))
                            warning(warnState);
                            % Check if we now have the correct function on
                            % the path
                            if ~strcmp(which(func2str(sceneGeometry.refraction.handle)), sceneGeometry.refraction.path)
                                % all set
                            else
                                error('Unable to re-instantiate the ray tracing function')
                            end
                        end
                    end
                end
                if verbose
                    fprintf('Compiled ray trace function available\n');
                end
            otherwise
                error('Unrecognized function type');
        end
    end
end

end