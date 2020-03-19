function glintCoord = obtainGlintCoord(imagePoints,pointLabels)
% Identify the image point that is the glint and return its coordinate
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
%
% Outputs:
%   glintCoord            - A nx2 vector with the image coordinates of the
%                           n glints.
%

% Find the pointLabels for which the first five characters are "glint".
% This handles finding the strings "glint_01", "glint_02", etc.
glintIdx = strncmp(pointLabels,'glint',5);

% If we have any such points, get the image coordinates and return
if any(glintIdx)
    glintCoord = imagePoints(glintIdx,:);
else
    glintCoord = [];
end

end % obtainGlintCoord

