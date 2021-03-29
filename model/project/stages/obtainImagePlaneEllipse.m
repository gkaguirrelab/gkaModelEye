function [pupilEllipse,pupilFitError,pointLabels] = obtainImagePlaneEllipse(imagePoints,pointLabels,sceneGeometry,p,eyePose)
% Fit an ellipse to the pupil image
%
% Syntax:
%  [pupilEllipse,pupilFitError,pointLabels] = obtainImagePlaneEllipse(imagePoints,pointLabels,sceneGeometry,p,eyePose)
%
% Description:
%   Fit an ellipse to the pupil in the image plane. This routine handles
%   both the simple case of a zero-thickness iris, and the more complex
%   case of an aperture stop that is modeled as having a thickness. In the
%   latter case, we must determine if the front or back edge of the iris
%   defines the pupil border.
%
% Inputs:
%   imagePoints           - nx2 vector. Points in image coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   p                     - Structure. The structure returned by the
%                           parameter parser in the calling function.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%
% Outputs:
%   pupilEllipse          - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
%   pupilFitError         - The RMSE distance (in pixels) of the pupil
%                           perimeter points to the pupil ellipse. This
%                           value indicates the degree to which the shape
%                           of the entrance pupil departs from an ellipse.
%   pointLabels           - nx1 cell array. The name of each eye point.
%

% If we have zero aperture stop radius, return
if ~(eyePose(4) > 0)
    pupilFitError = nan;
    pupilEllipse=nan(1,5);
    return
end

% Handle the simple case of a zero-thickness aperture stop
if sceneGeometry.eye.iris.thickness==0
    % The simple case of a zero-thickness aperture stop. Identify the
    % perimeter points.
    pupilPerimIdx = logical(strcmp(pointLabels,'pupilPerimeter'));
    % Fit the ellipse
    [pupilEllipse, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
    % We are done and can return
    return
end

% Handle the more complicated case of a non-zero thickness aperture stop.
% Identify the perimeter points for the front and back virtual image of the
% aperture stop.
pupilPerimIdxFront = logical(strcmp(pointLabels,'pupilPerimeterFront'));
pupilPerimIdxBack = logical(strcmp(pointLabels,'pupilPerimeterBack'));

% Obtain the ellipse fit to the front and back
frontEllipse = pupilEllipseFit(imagePoints(pupilPerimIdxFront,:));
backEllipse = pupilEllipseFit(imagePoints(pupilPerimIdxBack,:));

% If either ellipse fit yielded nans, then the final ellipse is nan
if any(isnan(frontEllipse)) || any(isnan(backEllipse))

    % If we are unable to fit an ellipse to the front or back pupil
    % perimeter, then exit with nans
    pupilFitError = nan;
    pupilEllipse=nan(1,5);

else
    
    % For each position on the perimeter of the pupil, determine which
    % point (front or back) is farther from the center of the ellipse.
    % These will be the hidden points
    centerDistance = sqrt(sum(((...
        bsxfun(@minus,...
        mean([frontEllipse(1:2);backEllipse(1:2)]),...
        imagePoints(logical(pupilPerimIdxFront+pupilPerimIdxBack),:))...
        ).^2),2));
    hideBack = (centerDistance(1:nStopPerimPoints)-centerDistance(nStopPerimPoints+1:nStopPerimPoints*2))<0;

    % Identify the set of non-hidden pupil perimeter points
    backStopVisibleIdx = pupilPerimIdxBack;
    backStopVisibleIdx(pupilPerimIdxBack)=~hideBack;
    frontStopVisibleIdx = pupilPerimIdxFront;
    frontStopVisibleIdx(pupilPerimIdxFront)=hideBack;
    pupilPerimIdx = or(backStopVisibleIdx,frontStopVisibleIdx);

    % Remove those pupil perimeter points that have had poor ray tracing
    pupilPerimIdx = and(pupilPerimIdx,targetIntersectError<p.Results.rayTraceErrorThreshold);

    % Fit the ellipse
    [pupilEllipse, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));

    % Update the pointLabels to indicate the hidden points
    pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)),'_hidden');
    pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)),'_hidden');

end


end % obtainImagePlaneEllipse