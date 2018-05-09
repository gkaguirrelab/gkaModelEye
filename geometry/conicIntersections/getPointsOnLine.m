% getPointsOnLine -  given an homogeneous line return two homogeneous points on it
%              - Pierluigi Taddei (pierluigi.taddei@polimi.it)
%
% Usage:    [p1 p2] = getPointsOnLine(l)
%
% Arguments:
%           l - line vector
%           p1 p2 - two points on l
%
% 06.3.2007 : Created
%
%
function [p1 p2] = getPointsOnLine(l)
   if (l(1) == 0 && l(2) == 0) %line at infinity
       p1 = [1 0 0]';
       p2 = [0 1 0]';
   else
    p2 = [-l(2), l(1), 0]';
    if (abs(l(1)) < abs(l(2)))
        p1 = [0, -l(3), l(2)]';
    else
        p1 = [-l(3), 0, l(1)]';
    end  
   end
end