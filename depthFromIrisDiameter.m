function cameraTranslationDepth = depthFromIrisDiameter( sceneGeometry, observedIrisDiamPixels )
% Estimate camera depth given sceneGeometry and iris diameter in pixels
%
% Syntax:
%  cameraTranslationDepth = depthFromIrisDiameter( sceneGeometry, observedIrisDiamPixels )
%
% Description:
%   There is limited individual variation in the horizontal visible
%   diameter of the human iris. The maximum observed diameter of the border of the
%   iris in a set of images of the eye will correspond to the diameter of
%   the iris when the eye is posed so a line that connects the center of
%   rotation of the eye and the center of the iris is aligned with the
%   optical axis of the camera. Given an observed iris diameter in pixels
%   in the image and a candidate sceneGeometry structure, we can calculate
%   the distance of the camera from the corneal surface of the eye.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry file for which we wish to
%                           refine the camera depth value.
%   observedIrisDiamPixels - Scalar. The maximum observed diameter of the
%                           iris across image frames, in units of pixels.
%
% Outputs;
%   cameraTranslationDepth - Scalar. The distance, in mm, of the camera
%                           from the corneal apex at pose angle [0, 0, 0].
%
% Examples:
%{
    %% Recover a veridical camera distance
    % Create a sceneGeometry structure
    sceneGeometry = createSceneGeometry();
    % Calculate what the observed iris diameter should be at 100 mm
    veridicalSceneGeometry = sceneGeometry;
    veridicalSceneGeometry.extrinsicTranslationVector(3) = 100;
    [~, imagePoints, ~, ~, pointLabels] = ...
    	pupilProjection_fwd([0 0 0 1], veridicalSceneGeometry, 'fullEyeModelFlag', true, 'nIrisPerimPoints', 100);
    idx = find(strcmp(pointLabels,'irisPerimeter'));
    observedIrisDiamPixels = max(imagePoints(idx,1))-min(imagePoints(idx,1));
    % Now call the estimation function with the default (incorrect) scene
    % geometry
    cameraTranslationDepth = ...
        depthFromIrisDiameter( sceneGeometry, observedIrisDiamPixels );
    % Report the results
    fprintf('Veridical camera depth: %4.0f, recovered camera depth: %4.0f \n',veridicalSceneGeometry.extrinsicTranslationVector(3),cameraTranslationDepth);
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('sceneGeometry', @isstruct);
p.addRequired('maxIrisDiameterPixels',@isnumeric);

% parse
p.parse(sceneGeometry, observedIrisDiamPixels)

x0 = sceneGeometry.extrinsicTranslationVector(3);

cameraTranslationDepth = fminsearch(@objfun, x0);
    function fVal = objfun(x)
        candidateSceneGeometry = sceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector(3) = x;
        [~, imagePoints, ~, ~, pointLabels] = ...
            pupilProjection_fwd([0 0 0 1], candidateSceneGeometry, 'fullEyeModelFlag', true, 'nIrisPerimPoints', 100);
        idx = find(strcmp(pointLabels,'irisPerimeter'));
        predictedIrisDiamPixels = max(imagePoints(idx,1))-min(imagePoints(idx,1));
        fVal = (predictedIrisDiamPixels - observedIrisDiamPixels)^2;
    end

end

