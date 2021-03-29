function classString = classify(S)
% Return a char vector that describes the form of the quadric surface
%
% Syntax:
%  classString = quadric.classify(S)
%
% Description:
%   Given a quadric surface (S) the routine returns a character vector that
%   identifies the quadric type.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%
% Outputs:
%   classString           - Character vector.
%

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% Place the quadric in canonical form
S = quadric.translate(S,-quadric.center(S));
S = quadric.alignAxes(S);

% obtain the eigenvalues
[~,evals] = svd(S( 1:3, 1:3 ) / -S( 4, 4 ));

% Identify the class of the quadric
switch num2str(sign(diag(evals))','%d')
    case '111'
        classString = 'ellipsoid';
    case {'11-1', '1-11', '-111'}
        classString = 'hyperboloid of one sheet';
    case {'1-1-1', '-11-1', '-1-11'}
        classString = 'hyperboloid of two sheets';
    otherwise
        classString = 'degenerate case';
end

end

