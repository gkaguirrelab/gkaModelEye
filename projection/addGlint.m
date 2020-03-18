function [eyePoints, pointLabels] = addGlint(eyePoints,pointLabels,sceneGeometry,p,eyePose)


%% Calc glint
% Find the location of the glint in the eye world coordinate frame. The
% glint is the reflection of a light source from the tear film of the eye.
% The location of the glint in the image is subject to refraction by
% artificial lenses.

%% Extract fields from Results struct
glintRayFunc = p.Results.glintRayFunc;
rayTraceErrorThreshold = p.Results.rayTraceErrorThreshold;

% Exit the routine if we do not have what we need
if ~isfield(sceneGeometry,'refraction') || isempty(glintRayFunc)
    return
end

if ~isfield(sceneGeometry.refraction,'glint')
    return
end




% The position of the light source in the world coordinate frame
% that is the source of the glint
glintSourceWorld = sceneGeometry.cameraPosition.translation + ...
    sceneGeometry.cameraPosition.glintSourceRelative;

% How many light sources do we have?
nGlints = size(glintSourceWorld,2);

% Assemble the args
args = {sceneGeometry.cameraPosition.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.refraction.cameraToMedium.opticalSystem, ...
    sceneGeometry.refraction.glint.opticalSystem, ...
    sceneGeometry.refraction.mediumToCamera.opticalSystem};

% Loop through the glints
for gg = 1:nGlints
    
    % Perform the computation using the passed function handle
    [virtualImageRay, ~, intersectError] = ...
        glintRayFunc(glintSourceWorld(:, gg), eyePose, args{:});
    
    % If we have a good trace, add the glint point and label
    if intersectError < rayTraceErrorThreshold
        eyePoints = [eyePoints; virtualImageRay(1,:)];
        glintLabel = sprintf('glint_%02d',gg);
        pointLabels = [pointLabels; {glintLabel}];
    end
    
end


end