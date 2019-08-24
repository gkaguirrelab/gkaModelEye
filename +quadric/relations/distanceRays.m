function [centerPoint, distance, R1p, R2p]=distanceRays(R1,R2)
% Finds the point that is mutually closest to two rays
%
% Syntax:
%  [centerPoint, distance, R1p, R2p] = quadric.distanceRays(R1,R2)
%
% Description:
%   Given two rays in standard form, returns the point that is mutually
%   closest to each of the rays. Adapted from code provided by Alexander
%   Brodsky on MATLAB central.
%
% Inputs:
%   R1, R2                - 3x2 matrix that specifies a vector of the form
%                           [p; d]:
%                               R = p + d,
%                           where p is vector origin, d is the direction.
%
% Outputs:
%   centerPoint           - Scalar. Point that is mutually closest to the
%                           two rays
%   distance              - Scalar. Distance between the two rays at their
%                           point of closest approach
%   R1p, R2p              - 3x1 arrays that provide the point on each ray
%                           that is closest to the other ray
%
% Examples:
%{
%}

P0=R1(:,1)';
P1=(R1(:,1)+R1(:,2))';
Q0=R2(:,1)';
Q1=(R2(:,1)+R2(:,2))';
u=P1-P0; v=Q1-Q0; w0=P0-Q0;
a=u*u'; b=u*v'; c=v*v'; d=u*w0'; e=v*w0';
sc=(b*e-c*d)/(a*c-b^2);
tc=(a*e-b*d)/(a*c-b^2);
distance=norm(w0+(sc*u-tc*v));
R1p=(P0+sc*u)';
R2p=(Q0+tc*v)';
centerPoint = (R1p+R2p)./2;

end