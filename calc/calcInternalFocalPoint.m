function [focalPoint,raySeparationAtFocalPoint,rayPath1,rayPath2] = calcInternalFocalPoint(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
% Focal point for rays arising from a specified field position
%
% Syntax:
%  [focalPoint,raySeparationAtFocalPoint,rayPath1,rayPath2] = calcInternalFocalPoint(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
%
% Description
%   This routine obtains the point of intersection of a pair of emergent
%   rays that are produced by a pair of diverging (or parallel) incident
%   rays. The incident rays arise from a particular field position (defined
%   w.r.t. the referenceCoord). If the rayOriginDistance is Inf, then the
%   rays are parallel. Otherwise, the rays diverge from a common origin
%   point at the specified distance. The point of closest approach of the
%   emergent rays is the internal focal point for the ray origin position.
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
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in referenceCoord.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the longitudinal axis origin.
%   referenceCoord        - 3x1 vector that provides the coordinate from
%                           which the ray origin angles and distance are
%                           to be calculated. By default, this is [0;0;0],
%                           which is the origin coordinate on the
%                           longitudinal axis.
%   rayIntersectionHeight - Scalar. The divergent rays will arrive at the
%                           corneal apex separated by 2x this value.
%   effectiveInfinity     - Scalar. Rays arising from this point or beyond
%                           will be treated as parallel.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   focalPoint            - 3x1 matrix with the Cartesian coordinates of
%                           point of closest approach of the two rays after
%                           emerging from the optical system. Ideally, this
%                           point is on the retina.
%   raySeparationAtFocalPoint - The distance between the two rays at their
%                           point of closest approach. Ideally, this value
%                           is zero.
%   rayPath1, rayPath2    - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%
% Examples:
%{
    eye = modelEyeParameters();
    [focalPoint, raySeparationAtFocalPoint] = calcInternalFocalPoint(eye);
%}

%% Arguments
arguments
    opticalSystem
    fieldAngularPosition (2,1) double = [0, 0]
    rayOriginDistance (1,1) double = Inf
    angleReferenceCoord (3,1) double = [0, 0, 0]
    distanceReferenceCoord (3,1) double = [0, 0, 0]
    rayIntersectionHeight (1,1) double = 0.25
    effectiveInfinity (1,1) double = 1e4
    cameraMedium = 'air'
end

% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium,...
            'opticalSystemNumRows',[]);
    end
end

% Define the rayOrigin. First, create a ray that leaves the
% angleReferenceCoord at the specified fieldAngularPositions
R = quadric.anglesToRay(angleReferenceCoord,fieldAngularPosition(1),fieldAngularPosition(2));

% We now need to determine how far we can travel along this ray such that
% our position from the distanceReferenceCoord is equal to the desired
% rayOriginDistance. I brute force this.
finiteDistance = min([rayOriginDistance effectiveInfinity]);
myObj = @(d) finiteDistance - norm( R(:,1)+R(:,2)*d-distanceReferenceCoord );
d = fzero(myObj,finiteDistance);

% Produce the specified rayOrigin
rayOrigin = R(:,1)+R(:,2)*d;

% The separation between the rays at the origin of the longitudinal axis.
deltaPosition = [0;rayIntersectionHeight/sqrt(2);rayIntersectionHeight/sqrt(2)];

% Create rays that start at rayOrigin and will intersect the xy plane
% of the longitudinal axis at a distance of 2 x rayIntersectionHeight.
if norm(rayOrigin) >= (effectiveInfinity-1e-6)
    
    % These rays are parallel
    myR1 = quadric.coordsToRay([rayOrigin+deltaPosition,angleReferenceCoord+deltaPosition]);
    myR2 = quadric.coordsToRay([rayOrigin-deltaPosition,angleReferenceCoord-deltaPosition]);
    
else
    
    % These rays are diverging
    myR1 = quadric.coordsToRay([rayOrigin,angleReferenceCoord+deltaPosition]);
    myR2 = quadric.coordsToRay([rayOrigin,angleReferenceCoord-deltaPosition]);
    
end

% Trace the rays
[outputRay1, rayPath1] = rayTraceQuadrics(myR1, opticalSystem);
[outputRay2, rayPath2] = rayTraceQuadrics(myR2, opticalSystem);

% The point of intersection of the rays within the eye
[focalPoint, raySeparationAtFocalPoint] = ...
    quadric.distanceRays(outputRay1,outputRay2);


end