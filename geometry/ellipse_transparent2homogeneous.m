function M = ellipse_transparent2homogeneous( t )
% Converts transparent ellipse parameters into homogeneous matrix form
%
% Syntax:
%  M = ellipse_transparent2homogeneous( t )
%
% Description:
%   The homogeneous representation of a conic is a matrix:
%
%       m = [A C D; C B E; D E F]
%
%   that represents the equation:
%
%       A x^2 + B y^2 + 2C xy + 2D x + 2Ey + F = 0 
%
%   The homogeneous form is used in routines that find the intersection
%   betweem conic sections (see: conicIntersections).
%
%   The transparent parameters are converted to explicit form (center x and
%   y, semi major and minor lengths, and tilt in radians) and then into
%   homogeneous form. The latter transformation makes use of equations
%   provided by Jan Benes.
%

% Convert from transparent to explicit form
A = GtoA(ellipse_transparent2ex(t));

% Convert from explicit to homogeneous form. Conversion formulae taken
% from: https://math.stackexchange.com/questions/159095/how-to-derive-ellipse-matrix-for-general-ellipse-in-homogenous-coordinates
% NOTE: The notation used on the linked web page switches the assignment of
% the B and C components compared to what is used here (and to what is used
% by Pierluigi Taddei in his conic intersection code).
%
% A x^2 + B y^2 + 2C xy + 2D x + 2Ey + F = 0
%


% Assemble into the homogeneous matrix
M = [A C D; C B E; D E F];

end

