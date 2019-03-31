function p = ellipse_ex2im(varargin)
% Cast ellipse defined with explicit parameters to implicit form.
% In explicit form, the parameters of the ellipse are its center (cx,cy) its semi-major
% and semi-minor axes (a and b) and its angle of tilt (theta).
%
% See also: ellipse_im2ex

% Copyright 2011 Levente Hunyadi

if nargin > 1
    narginchk(5,5);
    for k = 1 : 5
        validateattributes(varargin{k}, {'numeric'}, {'real','scalar'});
    end
    p = ellipse_implicit(varargin{:});
else
    narginchk(1,1);
    q = varargin{1};
    validateattributes(q, {'numeric'}, {'real','vector'});
    q = q(:);
    validateattributes(q, {'numeric'}, {'size',[5 1]});
    p = ellipse_implicit(q(1), q(2), q(3), q(4), q(5));
end

function p = ellipse_implicit(c1,c2,a,b,theta)
% Cast ellipse defined in Kepler form to standard parameter vector form.

p = imconictranslate(imconicrotate([1/a^2 0 1/b^2 0 0 -1], theta), [c1;c2]);
p = p ./ norm(p);