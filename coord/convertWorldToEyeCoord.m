function eyeCoord = convertWorldToEyeCoord(worldCoord)
% Convert a coordinate from eye space to world space
%
% Syntax:
%  eyeCoord = convertWorldToEyeCoord(worldCoord)
%
% Description
%   Separate coordinate systems are used for the "eye" and the "world".
%   This routine takes a world coordinate and returns the eye coordinate.
%
% Inputs:
%   worldCoord            - A 3x1 vector that specifies a point in world 
%                           coordinates (x, y, z).
%
% Outputs:
%   eyeCoord              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%

% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. We also transpose as eye coordinates are row vectors
% and world coordinates are column vectors
eyeCoord = worldCoord([3 1 2])';

end