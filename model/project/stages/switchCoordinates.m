function worldPoints = switchCoordinates(headPoints)
% Switch from head coordiates to world coordinates
%
% Syntax:
%  worldPoints = switchCoordinates(headPoints)
%
% Description:
%   After being subject to rotation, the eye points are in head coordinates
%   (which is essentially still the eye coordinate frame, just now defined
%   relative to the unrotated eye). We switch here to the world coordinate
%   frame. While the eye and world coordinates differ in their row/column
%   orientaton, we set the worldPoints matrix to match the orientation of
%   the headPoints and eyePoints.
%
% Inputs:
%   headPoints            - nx3 vector. Points in head coordinates.
%
% Outputs:
%   worldPoints           - nx3 vector. Points in world coordinates.
%

% Note the transpose operation to keep the worldPoints matrix in the same
% orientation as the source headPoints
worldPoints = convertEyeToWorldCoord(headPoints)';

end