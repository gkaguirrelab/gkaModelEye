function p = ellipsefit_robust(R, Q)
% An overloaded version of the ellipsefit_robust routine from quadfit
%
% Syntax:
%  p = ellipsefit_robust(R, Q)
%
% Description:
%   To support code generation, the computation in line 50 required
%   explicit conversion of the right hand side of the equation from complex
%   to real.
%
% The original header follows:


% Constrained ellipse fit by solving a modified eigenvalue problem.
% The method is numerically stable.
%
% Input arguments:
% R:
%    positive semi-definite data covariance matrix
% Q:
%    constraint matrix in parameters x^2, xy, y^2, x, y and 1.
%
% Output arguments:
% p:
%    estimated parameters (taking constraints into account)

% References:
% Radim Halir and Jan Flusser, "Numerically stable direct least squares fitting of
%    ellipses", 1998

% Copyright 2012 Levente Hunyadi

validateattributes(R, {'numeric'}, {'real','2d','size',[6,6]});
validateattributes(Q, {'numeric'}, {'real','2d','size',[6,6]});

% check that constraint matrix has all zeros except in upper left block
assert( nnz(Q(4:6,:)) == 0 );
assert( nnz(Q(:,4:6)) == 0 );

S1 = R(1:3,1:3);     % quadratic part of the scatter matrix
S2 = R(1:3,4:6);     % combined part of the scatter matrix
S3 = R(4:6,4:6);     % linear part of the scatter matrix
T = -(S3 \ S2');     % for getting a2 from a1
M = S1 + S2 * T;     % reduced scatter matrix
M = Q(1:3,1:3) \ M;  % premultiply by inv(C1), e.g. M = [M(3,:)./2 ; -M(2,:) ; M(1,:)./2] for an ellipse
[evec,~] = eig(M);   % solve eigensystem

% evaluate a'*C*a, e.g. cond = 4 * evec(1,:).*evec(3,:) - evec(2,:).^2 for an ellipse
cond = zeros(1,size(evec,2));
for k = 1 : numel(cond)
    cond(k) = real(evec(:,k)'*Q(1:3,1:3)*evec(:,k));
end

% eigenvector for minimum positive eigenvalue
evec = evec(:,cond > 0);
cond = cond(cond > 0);
[~,ix] = min(cond,[],1);
p1 = evec(:,ix);  % eigenvector for minimum positive eigenvalue

% ellipse coefficients
p = [p1 ; T * p1];
