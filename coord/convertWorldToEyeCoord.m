function eyeCoord = convertWorldToEyeCoord(worldCoord)
% Convert a coordinate from eye space to world space
%
% Syntax:
%  eyeCoord = convertWorldToEyeCoord(worldCoord)
%
% Description
%   Separate coordinate systems are used for the "eye" and the "world"


% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. We also transpose as eye coordinates are row vectors
% and world coordinates are column vectors
eyeCoord = worldCoord([3 1 2])';

end