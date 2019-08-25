function S = unitSphere
% Returns a quadric surface matrix that defines a unit sphere
%
% Syntax:
%  S = quadric.unitSphere
%
% Description:
%   Returns a quadric surface matrix that is a sphere with semi-radii of 1
%   unit.
%
% Inputs:
%   none
%
% Outputs:
%   S                     - 4x4 matrix of the quadric surface.
%

S = eye(4);
S(4,4)=-1;

end