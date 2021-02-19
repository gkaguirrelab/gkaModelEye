function [focalPoint, blur] = calcInternalFocalPoint(opticalSystem,rayOriginDistance,rayHeight)
% Focal point for rays arising at the specified point on the optical axis
%
% Syntax:
%  focalPoint = calcInternalFocalPoint(opticalSystem,rayOriginDistance,rayHeight,cameraMedium)
%
% Description
%   This routine examines the behavior of a pair of diverging rays that
%   arise on the optical axis at the specified rayOriginDistance. These
%   rays intersect the plane of the front surface of the optical system,
%   and ultimately give rise to a pair of emergent rays. The closest point
%   of approach of these emergent rays is the internal focal point for this
%   near point.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex. Assumed to be
%                           500 mm if not defined.
%
% Outputs:
%   focalPoint            - 3x1 matrix with the Cartesian coordinates of
%                           point of closest approach of the two rays after
%                           emerging from the optical system. Ideally, this
%                           point is on the retina.
%   blur                  - The distance between the two rays at their
%                           point of closest approach. Ideally, this value
%                           is zero.
%
% Examples:
%{
    eye = modelEyeParameters();
    opticalSystem = assembleOpticalSystem(eye,...
        'surfaceSetName','mediumToRetina','cameraMedium','air');
    [focalPoint, blur] = calcInternalFocalPoint(opticalSystem,Inf);
%}


% Handle missing inputs
if nargin<2
    error('calcInternalFocalPoint:invalidArguments','Too few input arguments')
end

if nargin==2
    rayHeight = 1;
end

% Create rays that start on the optical axis at the rayOriginDistance, and
% intersect the plane of the front surface of the cornea at the height
% given by ±rayHeight.
if isinf(rayOriginDistance)
    % The rays are fixed at parallel
    effectiveInfinity = 2000;
    myR1 = quadric.normalizeRay([effectiveInfinity,-1;rayHeight,0;0,0]);
    myR2 = quadric.normalizeRay([effectiveInfinity,-1;-rayHeight,0;0,0]);
else
    % The principal point of the optical system.
    myPrincipalPoint = calcPrincipalPoint(opticalSystem);
    
    % Account for the depth of the principal point in calculating the
    % position from which the rays arise, as the coordinate space is w.r.t.
    % the front corneal surface.
    myRayOrigin = (rayOriginDistance) - sum(myPrincipalPoint.*[1;0;0]);
    
    % Calculate the angle with which the rays diverge from the optical axis
    % such that they will intersect the principal plane at the ray height
    myAngle = rad2deg(atan2(rayHeight,-myRayOrigin));
    
    % Define the two rays
    myR1 = quadric.normalizeRay(quadric.anglesToRay([myRayOrigin;0;0],myAngle,0));
    myR2 = quadric.normalizeRay(quadric.anglesToRay([myRayOrigin;0;0],-myAngle,0));
end

% The point of intersection of the rays within the eye
[focalPoint, blur] = quadric.distanceRays( ...
    rayTraceQuadrics(myR1, opticalSystem), ...
    rayTraceQuadrics(myR2, opticalSystem)  );

end