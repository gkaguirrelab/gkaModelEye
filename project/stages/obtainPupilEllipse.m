function [pupilEllipseOnImagePlane,pupilFitError,pointLabels] = obtainPupilEllipse(imagePoints,pointLabels,sceneGeometry,p,eyePose)



%% Obtain the pupil ellipse
% Proceed with fitting if we have a non-zero stop radius
if eyePose(4) > 0
    if sceneGeometry.eye.iris.thickness==0
        % The simple case of a zero-thickness aperture stop. Identify the
        % perimeter points.
        pupilPerimIdx = logical(strcmp(pointLabels,'pupilPerimeter'));
        % Fit the ellipse
        [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
    else
        % The more complicated case of a non-zero thickness aperture stop.
        % Identify the perimeter points for the front and back virtual
        % image of the aperture stop.
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
            pupilEllipseOnImagePlane=nan(1,5);
        else
            % For each position on the perimeter of the pupil, determine
            % which point (front or back) is farther from the center of the
            % ellipse. These will be the hidden points
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
            % Remove those pupil perimeter points that have had poor ray
            % tracing
            pupilPerimIdx = and(pupilPerimIdx,targetIntersectError<p.Results.rayTraceErrorThreshold);
            % Fit the ellipse
            [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints(pupilPerimIdx,:));
            % Update the pointLabels to indicate the hidden points
            pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxBack,~backStopVisibleIdx)),'_hidden');
            pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)) = strcat(pointLabels(and(pupilPerimIdxFront,~frontStopVisibleIdx)),'_hidden');
        end
    end
else
    pupilFitError = nan;
    pupilEllipseOnImagePlane=nan(1,5);
end

end