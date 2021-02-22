function [focalPoint, raySeparationAtFocalPoint, rayPath1, rayPath2] = calcInternalFocalPoint(opticalSystem,rayOrigin,rayIntersectionHeight,effectiveInfinity)
% Focal point for rays arising at the specified point on the optical axis
%
% Syntax:
%  [focalPoint, raySeparationAtFocalPoint] = calcInternalFocalPoint(opticalSystem,rayOrigin,rayIntersectionHeight)
%
% Description
%   This routine examines the behavior of a pair of diverging rays that
%   arise from a specified location and are directed to either side of the
%   incident node of the optical system. The closest point of approach of
%   the resulting emergent rays is the internal focal point for this
%
% Inputs:
%   opticalSystem         - Either an eye structure (from which a
%                           "mediumToRetina" optical system in air will be
%                           derived), or an mx19 matrix, where m is set by
%                           the key value opticalSystemNumRows. Each row
%                           contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%   rayOrigin             - Scalar or 3x1 vector. If a scalar, this is the
%                           distance (in mm) from the corneal apex of a
%                           rayOrigin that is on the optical axis. A 3x1
%                           vector can specify a position that is not on
%                           the optical axis. A value of inf in the first
%                           or only element defines parallel rays arriving
%                           from effective infinity.
%   rayIntersectionHeight - Scalar. The divergent rays will arrive at the
%                           corneal apex separated by 2x this value.
%   effectiveInfinity     - Scalar. Rays arising from this point or beyond
%                           will be treated as parallel.
%
% Outputs:
%   focalPoint            - 3x1 matrix with the Cartesian coordinates of
%                           point of closest approach of the two rays after
%                           emerging from the optical system. Ideally, this
%                           point is on the retina.
%   raySeparationAtFocalPoint - The distance between the two rays at their
%                           point of closest approach. Ideally, this value
%                           is zero.
%
% Examples:
%{
    eye = modelEyeParameters();
    [focalPoint, raySeparationAtFocalPoint] = calcInternalFocalPoint(eye,Inf);
%}


% Handle missing inputs
if nargin<2
    error('calcInternalFocalPoint:invalidArguments','Too few input arguments')
end

if nargin==2
    rayIntersectionHeight = 0.5;
    effectiveInfinity = 2000;
end

if nargin==3
    effectiveInfinity = 2000;
end

% If rayOrigin is a scalar, place it on the optical axis
if isscalar(rayOrigin)
    if isinf(rayOrigin)
        rayOrigin = [effectiveInfinity;0;0];
    else
        rayOrigin = [rayOrigin;0;0];
    end
else
    if any(isinf(rayOrigin))
        error('calcInternalFocalPoint:invalidArguments','Provide a finite rayOrigin')
    end
end

% Force rayOrigin to be a column vector
if all(size(rayOrigin)==[1 3])
    rayOrigin = rayOrigin';
end

% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium','air',...
            'opticalSystemNumRows',[]);
    end
end

% Create rays that start at rayOrigin and diverge such that they will be
% separated by 2 x rayIntersectionHeight when they arrive at the plane of
% the corneal apex.
if norm(rayOrigin) >= (effectiveInfinity-1e-6)
        
    % Find the angle at which these rays will be traveling towards the
    % origin
    [p1p2, p1p3] = quadric.rayToAngles( quadric.normalizeRay([rayOrigin,-rayOrigin]) );
    
    % Define a delta separating of the rays
    deltaPosition = [0;rayIntersectionHeight;rayIntersectionHeight];
    
    % The rays are fixed at parallel
    myR1 = quadric.normalizeRay(quadric.anglesToRay(rayOrigin+deltaPosition,p1p2,p1p3));
    myR2 = quadric.normalizeRay(quadric.anglesToRay(rayOrigin-deltaPosition,p1p2,p1p3));
    
else
    
    % Find the angle at which these rays will be traveling towards the
    % origin
    [p1p2, p1p3] = quadric.rayToAngles( quadric.normalizeRay([rayOrigin,-rayOrigin]) );
    
    % Calculate the angle with which the rays must diverge at the rayOrigin
    % such that they are separated by 2 x rayIntersectionHeight in each
    % plane (p1p2 and p1p3) when they reach origin of the coordinate system
    deltaAngle = rad2deg(atan2(rayIntersectionHeight,norm(rayOrigin)));
    
    % Define the two rays.
    myR1 = quadric.normalizeRay(quadric.anglesToRay(rayOrigin,p1p2+deltaAngle,p1p3+deltaAngle));
    myR2 = quadric.normalizeRay(quadric.anglesToRay(rayOrigin,p1p2-deltaAngle,p1p3-deltaAngle));
    
end

% Trace the rays
[outputRay1, rayPath1] = rayTraceQuadrics(myR1, opticalSystem);
[outputRay2, rayPath2] = rayTraceQuadrics(myR2, opticalSystem);

% The point of intersection of the rays within the eye
[focalPoint, raySeparationAtFocalPoint] = ...
    quadric.distanceRays(outputRay1,outputRay2);

% Add the final step to the rayPath
rayPath1(:,end+1)=outputRay1(:,1);
rayPath2(:,end+1)=outputRay2(:,1);

end