function outP = ellipsefit_direct(x,y,ellipseFitRobustFunc)
% Direct least squares fitting of ellipses.
%
% Input arguments:
% x,y;
%    x and y coodinates of 2D points
%
% Output arguments:
% p:
%    a 6-parameter vector of the algebraic ellipse fit with
%    p(1)*x^2 + p(2)*x*y + p(3)*y^2 + p(4)*x + p(5)*y + p(6) = 0
%
% References:
% Andrew W. Fitzgibbon, Maurizio Pilu and Robert B. Fisher, "Direct Least
%    Squares Fitting of Ellipses", IEEE Trans. PAMI 21, 1999, pp476-480.

% Copyright 2011 Levente Hunyadi

validateattributes(x, {'numeric'}, {'real','nonempty','vector'});
validateattributes(y, {'numeric'}, {'real','nonempty','vector'});
x = x(:);
y = y(:);

% Determine if we have the compiled robust function available
if nargin == 2
    if exist('ellipseFit_robustMex','file') == 3
        ellipseFitRobustFunc = @ellipsefit_robustMex;
    else
        ellipseFitRobustFunc = @ellipsefit_robust;
    end
end

% normalize data
mx = mean(x);
my = mean(y);
sx = (max(x)-min(x))/2;
sy = (max(y)-min(y))/2;
smax = max(sx,sy);
sx = smax;
sy = smax;
x = (x-mx)/sx;
y = (y-my)/sy;

% build design matrix
D = [ x.^2  x.*y  y.^2  x  y  ones(size(x)) ];

% build scatter matrix
S = D'*D;

% build 6x6 constraint matrix
C = zeros(6,6);
C(1,3) = -2;
C(2,2) = 1;
C(3,1) = -2;

% fit
p = ellipseFitRobustFunc(S,-C);


% unnormalize
outP = zeros(6,1);
outP(1) = real(p(1)*sy*sy);
outP(2) = real(p(2)*sx*sy);
outP(3) = real(p(3)*sx*sx);
outP(4) = real(-2*p(1)*sy*sy*mx - p(2)*sx*sy*my + p(4)*sx*sy*sy);
outP(5) = real(-p(2)*sx*sy*mx - 2*p(3)*sx*sx*my + p(5)*sx*sx*sy);
outP(6) = real(p(1)*sy*sy*mx*mx + p(2)*sx*sy*mx*my + p(3)*sx*sx*my*my - p(4)*sx*sy*sy*mx - p(5)*sx*sx*sy*my + p(6)*sx*sx*sy*sy);

outP = outP ./ norm(outP);