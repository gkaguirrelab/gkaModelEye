function [Wp,rss] = quad2dproj(W, varargin)
% Projects a set of points onto a quadratic curve.
% The algorithm is proven to converge and reaches an accuracy of 7-8
% significant digits, and it takes 4-5 iterations per point, on average.
%
% Input arguments:
% W:
%    data points to project to the quadratic curve
%
% Output arguments:
% Wp:
%    data points that lie on the quadratic curve
% rss:
%    residual sum of squares (i.e. the sum of squares of the distances)

% References:
% Nikolai Chernov and H. Ma, "Least squares fitting of quadratic curves
%    and surfaces", Computer Vision, 2011, pages 285--302
% D. Eberly, "Distance from a point to an ellipse in 2D" (2004), Geometric
%    Tools, LLC, http://www.geometrictools.com

% Copyright 2013 Levente Hunyadi

% Currently implements projection only to ellipses

narginchk(2,4);
validateattributes(W, {'numeric'}, {'real','finite','2d'}, mfilename, 'W', 1);
if size(W,1) > size(W,2)
    validateattributes(W, {'numeric'}, {'size',[NaN,2]}, 1);
else
    validateattributes(W, {'numeric'}, {'size',[2,NaN]}, 1);
end
if nargin > 2
    narginchk(4,4);
    [center, axes, angle] = varargin{:};
    validateattributes(center, {'numeric'}, {'real','finite','vector','numel',2}, mfilename, 'center', 2);
    validateattributes(axes, {'numeric'}, {'real','finite','vector','numel',2}, mfilename, 'axes', 3);
    validateattributes(angle, {'numeric'}, {'real','finite','scalar'}, mfilename, 'angle', 4);
    center = center(:);
    axes = axes(:);
else
    narginchk(2,2);
    params = varargin{:};
    validateattributes(params, {'numeric'}, {'real','finite','vector','numel',5}, mfilename, 'p', 2);
    params = params(:);
    center = params(1:2);
    axes = params(3:4);
    angle = params(5);
end

if nargout > 1
    [Wp,rss] = ellipseproj(W, center, axes, angle);
else
    Wp = ellipseproj(W, center, axes, angle);
end