function objectHandle = addTransparentEllipseToFigure(pupilEllipseParams,imageSizeX,imageSizeY,Color,LineWidth,fImplicitPresent)
% Adds an ellipse in transparent params to the current figure
%
% Syntax:
%  objectHandle = addTransparentEllipseToFigure(pupilEllipseParams,imageSizeX,imageSizeY,Color,LineWidth,fImplicitPresent)
%
% Description:
%  Plots an ellipse in the currently active window and returns the handle
%  to the ellipse object
%
% Inputs:
%   pupilEllipseParams    - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
%   imageSizeX,imageSizeY - Scalar. Defines the image dimensions.
%   Color                 - Char vector. The color of the ellipse.
%   LineWidth             - Scalar. The line width.
%   fImplicitPresent      - Logical. Is the fimplicit function available?
%
% Outputs:
%   objectHandle          - Handle. The handle to the plotted ellipse.
%
% Examples:
%{
    pupilEllipseParams = [200, 200, 1000, 0.75, pi/4];
    figure
    addTransparentEllipseToFigure(pupilEllipseParams);
%}

% Handle incomplete inputs
switch nargin
    case 1
        imageSizeX = 400;
        imageSizeY = 400;
        Color = 'green';
        LineWidth = 1;
        fImplicitPresent = (exist('fimplicit','file')==2);
    case 2
        imageSizeY = imageSizeX;
        Color = 'green';
        LineWidth = 1;
        fImplicitPresent = (exist('fimplicit','file')==2);
    case 3
        Color = 'green';
        LineWidth = 1;
        fImplicitPresent = (exist('fimplicit','file')==2);
    case 4
        LineWidth = 1;
        fImplicitPresent = (exist('fimplicit','file')==2);
    case 5
        fImplicitPresent = (exist('fimplicit','file')==2);
    case 6
        % We have everything
    otherwise
        error('Not the right number of arguments for this function')
end

% Convert the ellipse parameters from transparent to implicit
pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilEllipseParams));

% Define the implicit function for the ellipse
fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);

% Superimpose the ellipse using fimplicit or ezplot (ezplot is the
% fallback option for older Matlab versions)
if fImplicitPresent
    objectHandle = fimplicit(fh,[1, imageSizeX, 1, imageSizeY],'Color', Color,'LineWidth',LineWidth);
else
    objectHandle = ezplot(fh,[1, imageSizeX, 1, imageSizeY]);
    set(objectHandle, 'Color', Color)
    set(objectHandle,'LineWidth',LineWidth);
end

end % function
