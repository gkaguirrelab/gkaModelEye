function glintCoord = obtainGlintCoord(imagePoints,pointLabels)
% Identify the image point that is the glint and return its coordinate
%
% Syntax:
%  glintCoord = obtainGlintCoord(imagePoints,pointLabels)
%
% Description:
%   This little routine finds the imagePoints that have a label that
%   contains the char vector "glint", and returns their coordinates to the
%   main function.
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

