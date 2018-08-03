function S = unitTwoSheetHyperboloid
% Returns a quadric surface matrix that is a unit two sheeted hyperboloid
%
% Syntax:
%  S = quadric.unitTwoSheetHyperboloid
%
% Description:
%   Returns a quadric surface matrix that is a two sheeted hyperboloid with
%   semi-radii of 1 unit. The surface opens along the first (x) dimension.
%
% Inputs:
%   none
%
% Outputs:
%   S                     - 4x4 matrix of the quadric surface.
%

S = eye(4);
S(1,1)=-1;
S(4,4)=1;

end