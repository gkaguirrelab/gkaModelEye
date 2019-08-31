function p = ellipse_ex2transparent(varargin)
% Cast ellipse defined in explicit form to transparent form
%
% Syntax:
%  p = ellipse_ex2transparent(explicitEllipseParams)
%
% Description:
%   In explicit form, the parameters of the ellipse are its center (cx,cy)
%   its semi-major and semi-minor axes (a and b) and its angle of tilt
%   (theta).
%
%   In transparent form, the parameters of the ellipse are its center
%   (cx,cy), its area (a), its eccentricity (e), and its angle of tilt
%   (theta).
%
% Inputs:
%   Either a single vector (5x1) of ellipse parameters in explicit form, or
%   five separate scalar variables, containing the five parameters.
%
% Outputs:
%   p                     - The parameters of an ellipse in transparent 
%                           form as a single, 5x1 vector
%
% Examples:
%{
    p1 = [100 100 400 0.9 pi/2];
    p2 = ellipse_ex2transparent(ellipse_transparent2ex(p1));
    assert(max(abs(p1-p2))<1e-6);
%}

if nargin > 1
    narginchk(5,5);
    for k = 1 : 5
        validateattributes(varargin{k}, {'numeric'}, {'real','scalar'});
    end
    p = ellipse_transparent(varargin{:});
else
    narginchk(1,1);
    q = varargin{1};
    validateattributes(q, {'numeric'}, {'real','vector'});
    q = q(:);
    validateattributes(q, {'numeric'}, {'size',[5 1]});
    p = ellipse_transparent(q(1), q(2), q(3), q(4), q(5));
end

end % function

%% Local function

function p = ellipse_transparent(c1,c2,a,b,theta)
% Cast ellipse defined in explict form to transparent parameter vector form.

p = nan(1,5);
p(1)=c1;
p(2)=c2;
p(3)=pi*a*b;
p(4)=sqrt(1- (b^2/a^2));
p(5)=theta;
end