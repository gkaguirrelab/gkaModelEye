function pt = imconictranslate(p, t)
% Translates a conic section given as an implicit equation.
%
% Input arguments:
% p:
%    the parameter vector p = [a b c d f g]
% t:
%    the translation vector t = [x y]
%
% See also: imconictranslation, imconicrotate

% Copyright 2010 Levente Hunyadi

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
assert(length(p) == 6, 'imconictranslate:DimensionMismatch', ...
    'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
if nargin > 1 && ~isempty(t)
    validateattributes(t, {'numeric'}, {'nonempty','real','vector'});
    assert(length(t) == 2, 'imconictranslate:DimensionMismatch', ...
        'A translation vector of length 2 expected.');
else
    t = imconictranslation(p);  % translate conic to origin
end

a = p(1); b = p(2); c = p(3); d = p(4); f = p(5); g = p(6);
x = t(1); y = t(2);

dt = d - 2*a*x - b*y;
ft = f - b*x - 2*c*y;
gt = g + a*x^2 + b*x*y + c*y^2 - d*x - f*y;

pt = zeros(size(p));
pt(:) = [a b c dt ft gt];
