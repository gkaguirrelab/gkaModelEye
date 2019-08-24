function S = unitParaboloid
% Returns a quadric surface matrix that is a unit paraboloid
%
% Syntax:
%  S = quadric.unitParaboloid
%
% Description:
%   Returns a quadric surface matrix that is a paraboloid with semi-radii
%   of 1 unit. The surface opens along the first (x) dimension.
%
% Inputs:
%   none
%
% Outputs:
%   S                     - 4x4 matrix of the quadric surface.
%

S = eye(4);
S(1,1)=0;
S(1,4)=-1;
S(4,1)=-1;
S(4,4)=0;

end