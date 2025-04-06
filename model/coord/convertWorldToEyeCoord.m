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
%	The eye coordinate frame is in mm and has the dimensions (p1,p2,p3).
%   The diagram is of a cartoon pupil, viewed directly from the front
%   (i.e., from the perspective of a camera looking at the eye of a
%   subject).
%
%                 |
%     ^         __|__
%  +  |        /     \
% p3  -  -----(   +   )-----
%  -  |        \_____/
%     v           |
%                 |
%
%           - <--p2--> +
%
%   Coordinate [0,0,0] corresponds to the apex (front surface) of the
%   cornea (when the apex of the cornea is aligned with the longitudinal
%   axis). The first dimension is depth, and has a negative value toward
%   the back of the eye. For the right eye, negative values on the p2
%   dimension are more temporal, and positive values are more nasal.
%   Positive values of p3 are upward, and negative values are downward.
%
% Inputs:
%   worldCoord            - 3xn vector that specifies points in world 
%                           coordinates (x, y, z).
%
% Outputs:
%   eyeCoord              - nx3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%

% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. We also transpose as eye coordinates are row vectors
% and world coordinates are column vectors
eyeCoord = worldCoord([3 1 2],:)';

end