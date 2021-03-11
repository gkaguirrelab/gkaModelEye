function [opticalSystem,surfaceLabels,surfaceColors] = parseOpticalSystemArgument(a,surfaceSet,cameraMedium)
% Returns an optical system for passed optical system, eye, or scene
%
% Syntax:
%  [opticalSystem,surfaceLabels,surfaceColors] = parseOpticalSystemArgument(a,surfaceSet,cameraMedium)
%
% Description:
%   For use in parsing the input arguments to a calling function
%

arguments
    a {mustBeOpticalSystem}
    surfaceSet char = 'mediumToRetina'
    cameraMedium = 'air'
end


% If the passed variable is a matrix, then just return it
if isnumeric(a)
    opticalSystem = a;
    surfaceLabels = {};
    surfaceColors = {};
    return
end

% If we were passed a structure, determine if it is a sceneGeometry or eye
% structure.
if isstruct(a)
    
    % Handle a surfaceSet
    if isfield(a,'opticalSystem')
        opticalSystem = a.opticalSystem;
        surfaceLabels = a.surfaceLabels;
        surfaceColors = a.surfaceColors;
        return
    end
    
    % Handle an eye structure
    if isfield(a,'cornea')
        [opticalSystem, surfaceLabels, surfaceColors] = ...
            assembleOpticalSystem(a,'surfaceSetName',surfaceSet,'cameraMedium',cameraMedium);
        return
    end
    
    % We have a sceneGeometry structure
    if isfield(a,'eye')
        % If we have a refraction field, see if we have this surface set
        if isfield(a,'refraction')
            if isfield(a.refraction,surfaceSet)
                opticalSystem = a.refraction.(surfaceSet).opticalSystem;
                surfaceLabels = a.refraction.(surfaceSet).surfaceLabels;
                surfaceColors = a.refraction.(surfaceSet).surfaceColors;
                return
            end
        end
        % Assemble the surface set from the eye field of sceneGeometry
        [opticalSystem, surfaceLabels, surfaceColors] = ...
            assembleOpticalSystem(a.eye,'surfaceSetName',surfaceSet,'cameraMedium',cameraMedium);
        return
    end
    
end

error('parseOpticalSystemArgument:noValidInterpretation','Could not create a valid opticalSystem from this input')

end