function [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints)
% Fit an ellipse to a set of image points
%
% Syntax:
%  [pupilEllipseOnImagePlane, pupilFitError] = pupilEllipseFit(imagePoints)
%
% Description:
%   This is a wrapper for the function ellipsefit_direct. The primary
%   application of this function is to fit an ellipse to a set of points in
%   an image that are on the boundary of the pupil. The purpose of this
%   function is to sanity check input, handle errors and warnings that can
%   arise during fitting, and provide a measure of fit error.
%
% Inputs:
%   imagePoints           - An nx2 matrix which contains the x and y
%                           positions of the n pupil perimeter points. At
%                           least five points are required.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form
%   pupilFitError         - The RMSE distance (in pixels) of the pupil
%                           perimeter points to the pupil ellipse.
%

% Make the codegen happy
coder.extrinsic('warning');

% Set up a variable to hold the pupil fit error
pupilFitError = nan;

% Before we try to fit the ellipse, make sure that the there are at least 5
% perimeter points that are non nan.
validPerimIdx = find(~any(isnan(imagePoints)')');

if ~isreal(imagePoints) || length(validPerimIdx)<5
    pupilEllipseOnImagePlane=nan(1,5);
else
    % Silence a warning that can arise regarding a nearly singular matrix
    
    warnState = warning;
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:singularMatrix');
    
    % We place the ellipse fit in a try-catch block, as the fit can fail
    % when the ellipse is so eccentric that it approaches a line. Try-catch
    % is not supported in compiled code, so we detect this situation and
    % proceed without the safety net.
    if isempty(coder.target)
        try
            % Ellipse fitting with routine from the quadfit toolbox
            implicitEllipseParams = real(ellipsefit_direct( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2)));
            % Convert the ellipse from implicit to transparent form
            pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
            % Place theta within the range of 0 to pi
            if pupilEllipseOnImagePlane(5) < 0
                pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
            end
            % Get the error of the ellipse fit to the pupil points
            pupilFitError = sqrt(nanmean(ellipsefit_distance( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2),ellipse_transparent2ex(pupilEllipseOnImagePlane)).^2));
        catch ME
            % If the ellipse fit direct fails because of an inability to
            % fit an ellipse to the provided points, return nans and issue
            % a warning.
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
    else
        % Ellipse fitting with routine from the quadfit toolbox
        implicitEllipseParams = real(ellipsefit_direct( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2)));
        % Convert the ellipse from implicit to transparent form
        pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
        % Place theta within the range of 0 to pi
        if pupilEllipseOnImagePlane(5) < 0
            pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
        end
        % Get the error of the ellipse fit to the pupil points
        pupilFitError = sqrt(nanmean(ellipsefit_distance( imagePoints(validPerimIdx,1), imagePoints(validPerimIdx,2),ellipse_transparent2ex(pupilEllipseOnImagePlane)).^2));
    end
    warning(warnState);
end
end % pupilEllipseFit