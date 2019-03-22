function [centerError, shapeError, areaError] = csaEllipseError(targetEllipse,fittedEllipse)
% Returns the mismatch in ellipse centers, shape, and area
%
% Syntax:
%  [centerError, shapeError, areaError] = csaEllipseError(targetEllipse,fittedEllipse)
%
% Description
%   Given a target and a fitted ellipse, this routine calculates three
%   error metrics that reflect differences between the ellipses.
%
% Inputs:
%   targetEllipse, fittedEllipse - 1x5 vectors that specify ellipses in
%                           transparent parameter form.
%
% Outputs:
%   centerError           - Scalar. The Euclidean distance between the
%                           ellipse centers.
%   shapeError            - Scalar. A value between 0 and unity that
%                           expresses the difference between the ellipses
%                           in their eccentricity and theta (tilt).
%   areaError             - Scalar. A value between 0 an unity that
%                           expresses the relative difference in area
%                           between the ellipses.
%

% The Euclidean distance between the ellipse centers
centerError = sqrt((targetEllipse(1) - fittedEllipse(1))^2 + ...
    (targetEllipse(2) - fittedEllipse(2))^2);

% The theta and eccentricity of an ellipse can be described as
% a point in polar coordinates. We express the error as
% the vector distance between these points. Direct minimization
% of differences in theta is a poor constraint, as differences
% in theta have reduced meaning at small eccentricities.
% Because ellipses are symmetric, theta spans the range of
% 0:pi. Therefore, the theta value is doubled prior to
% conversion to Cartesian coordinates so that the space wraps
% at the 0 - pi transition point. Eccentricity has a value
% ranging from zero (circular) to 1 (a fully flattened
% ellipse). I linearize the eccentricity value so that the
% error metric is sensitive to small differences in
% eccentricity. Finally, the value is divided by 2, so that the
% largest possible error is unity.
thetaT = targetEllipse(5)*2;
thetaC = fittedEllipse(5)*2;
rhoT = 1-sqrt(1-targetEllipse(4)^2);
rhoC = 1-sqrt(1-fittedEllipse(4)^2);
shapeError = ...
    sqrt(rhoT^2 + rhoC^2 - 2*rhoT*rhoC*cos(thetaT-thetaC))/2;

% Proportional difference in ellipse areas. Error value ranges from 0 to
% unity.
areaError = ...
    abs(targetEllipse(3) - fittedEllipse(3))/targetEllipse(3);

            
end

