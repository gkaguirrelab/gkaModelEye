function worldPoints = convertEyeToWorldCoord(eyePoints)
%% Switch axes to world coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (taken from Leszek Swirski), being
% viewed from above:
%
%    ^
%    |
%    |    .-.
% -Z |   |   | <- Head
%    +   `^u^'
% +Z |
%    |
%    |      W <- Camera    (As seen from above)
%    V
%
%     <-----+----->
%        -X   +X
%
% +X = right
% +Y = up
% +Z = front (towards the camera)
%
% The origin [0,0,0] corresponds to the front surface of the eye and the
% pupil center when the optical axis of the eye and the camera axis are
% aligned.

% Re-arrange the headPoints to transform to the world coordinate frame
worldPoints = eyePoints(:,[2 3 1]);

end