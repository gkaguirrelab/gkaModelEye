function varargout = ellipse_im2ex(varargin)
% Cast ellipse defined with implicit parameter vector to explicit form.
%
% See also: ellipse_ex2im

% Copyright 2011 Levente Hunyadi

if nargin > 1
    narginchk(6,6);
    for k = 1 : 6
        validateattributes(varargin{k}, {'numeric'}, {'real','scalar'});
    end
    [c1,c2,a,b,phi] = ellipse_explicit(varargin{:});
else
    narginchk(1,1);
    p = varargin{1};
    validateattributes(p, {'numeric'}, {'real','vector'});
    p = p(:);
    validateattributes(p, {'numeric'}, {'size',[6 1]});
    [c1,c2,a,b,phi] = ellipse_explicit(p(1), 0.5*p(2), p(3), 0.5*p(4), 0.5*p(5), p(6));
end
if nargout > 1
    varargout = num2cell([c1,c2,a,b,phi]);
else
    varargout{1} = [c1,c2,a,b,phi];
end

function [c1,c2,semia,semib,phi] = ellipse_explicit(a,b,c,d,f,g)
% Cast ellipse defined with explicit parameter vector to implicit form.

% helper quantities
N = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g);
D = b^2-a*c;
S = realsqrt((a-c)^2+4*b^2);

% semi-axes
ap = realsqrt( N/(D*(S-(a+c))) );
bp = realsqrt( N/(D*(-S-(a+c))) );
semia = max(ap,bp);
semib = min(ap,bp);

% center
c1 = (c*d-b*f)/D;
c2 = (a*f-b*d)/D;

% angle of tilt
if b ~= 0
    if abs(a) < abs(c)
        phi = 0.5*acot((a-c)/(2*b));
    else
        phi = 0.5*pi+0.5*acot((a-c)/(2*b));
    end
else
    if abs(a) < abs(c)
        phi = 0;
    else  % a > c
        phi = 0.5*pi;
    end
end