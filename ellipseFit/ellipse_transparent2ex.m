function p = ellipse_transparent2ex(varargin)
% Cast ellipse defined in transparent form to explicit form
%
% Syntax:
%  p = ellipse_transparent2ex(transparentEllipseParams)
%
% Description:
%   In transparent form, the parameters of the ellipse are its center
%   (cx,cy), its area (a), its eccentricity (e), and its angle of tilt
%   (theta).
%
%   In explicit form, the parameters of the ellipse are its center (cx,cy)
%   its semi-major and semi-minor axes (a and b) and its angle of tilt
%   (theta).
%
% Inputs:
%   Either a single vector (5x1) of ellipse parameters in transparent form,
%   or five separate scalar variables, containing the five parameters.
%
% Outputs:
%   p                     - The parameters of an ellipse in explicit form 
%                           as a single, 5x1 vector
%


if nargin > 1
    narginchk(5,5);
    for k = 1 : 5
        validateattributes(varargin{k}, {'numeric'}, {'real','scalar'});
    end
    p = ellipse_explicit(varargin{:});
else
    narginchk(1,1);
    q = varargin{1};
    validateattributes(q, {'numeric'}, {'real','vector'});
    q = q(:);
    validateattributes(q, {'numeric'}, {'size',[5 1]});
    p = ellipse_explicit(q(1), q(2), q(3), q(4), q(5));
end

end % function


%% LOCAL FUNCTION

function p = ellipse_explicit(c1,c2,area,e,theta)
% Cast ellipse defined in transparent form to explicit parameter vector form.

p(1)=c1;
p(2)=c2;
p(3)= sqrt(area / (pi * sqrt(1-e^2)));
p(4)= area / (pi*p(3));
p(5)=theta;
end