
function [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints)

% Set up a variable to hold the pupil fit error
pupilFitError = nan;

% Before we try to fit the ellipse, make sure that the radius is not zero,
% and that there are at least 5 perimeter points that are non nan.
validPerimIdx = find(~any(isnan(imagePoints)')');

if ~isreal(imagePoints) || length(validPerimIdx)<5
    pupilEllipseOnImagePlane=nan(1,5);
else
    % Silence a warning that can arise regarding a nearly singular matrix
    warnState = warning;
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:singularMatrix');

    % We place the ellipse fit in a try-catch block, as the fit can fail
    % when the ellipse is so eccentric that it approaches a line
    try
        % Ellipse fitting with routine from the quadfit toolbox
        implicitEllipseParams = ellipsefit_direct( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2));
        % Convert the ellipse from implicit to transparent form
        pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
        % Place theta within the range of 0 to pi
        if pupilEllipseOnImagePlane(5) < 0
            pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
        end
        % Get the error of the ellipse fit to the pupil points
        pupilFitError = sqrt(nanmean(ellipsefit_distance( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2),ellipse_transparent2ex(pupilEllipseOnImagePlane)).^2));
    catch ME
        % If the ellipse fit direct fails because of an inability to fit an
        % ellipse to the provided points, return nans and issue a warning.
        pupilEllipseOnImagePlane=nan(1,5);
        pupilFitError = nan;
        switch ME.identifier
            case {'MATLAB:badsubscript','MATLAB:realsqrt:complexResult','MATLAB:expectedReal','MATLAB:quad2dproj:expectedFinite','MATLAB:eig:matrixWithNaNInf'}
                warning('pupilProjection_fwd:ellipseFitFailed','Could not fit a valid ellipse to the pupil points; returning nans.');
            otherwise
                ME.identifier
                warning('pupilProjection_fwd:ellipseFitUnknownError','Undefined error during ellipse fitting to pupil perimeter; returning nans.');
        end
    end % try-catch block
    warning(warnState);
end
end % pupilEllipseFit