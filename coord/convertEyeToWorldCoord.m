function worldCoord = convertEyeToWorldCoord(eyeCoord)
% Convert a coordinate from world space to eye space
%
% Syntax:
%  eyeCoord = convertWorldToEyeCoord(worldCoord)
%
% Description
%   Separate coordinate systems are used for the "eye" and the "world".
%   This routine takes a world coordinate and returns the eye coordinate.
%
% Inputs:
%   eyeCoord              - nx3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%
% Outputs:
%   worldCoord            - 3xn vector that specifies points in world 
%                           coordinates (x, y, z).
%
% Examples:
%{
    % Confirm invertibility of the transform
    eyeCoord = rand(10,3);
    eyeCoordRecovered = convertWorldToEyeCoord(convertEyeToWorldCoord(eyeCoord));
    assert(max(max(eyeCoord-eyeCoordRecovered)) < eps);
%}

% Rearrange the eyeWorld dimensions to switch from eye to world
% coordinate space. We also transpose as eye coordinates are row vectors
% and world coordinates are column vectors
worldCoord = eyeCoord(:,[2 3 1])';

end