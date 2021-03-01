function [opticalPower, focalPoint] = calcOpticalPower(opticalSystem, rayOriginDistance, rayIntersectionHeight, cameraMedium)
% Calcuate the refractive power of an opticalSystem
%
% Syntax:
%  [opticalPower, focalPoint] = calcOpticalPower(opticalSystem, rayOriginDistance, rayIntersectionHeight, cameraMedium)
%
% Description
%   Calculates the refractive power of an optical system in units of
%	diopters. A negative value specifies the power of a lens that would be
%	worn by someone with myopia to correct their vision.
%
%   Some optical systems end in a medium with a refractive index other than
%   unity. In this case the optical power is given by:
%
%       diopters = refractiveIndex / effectiveFocalLength
%
%   where the effective focal length is the distance between the principal
%   point of the optical system and the focal point.
%
%   The calculation is performed along the longitudinal axis. Note that for
%   optical systems with spherical aberration, the calculated optical power
%   will vary depending upon the path of the ray. The rayIntersectionHeight
%   parameter controls the height at which the ray strikes the first
%   surface.
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
%   rayOriginDistance     - Scalar. Point of origin of the ray along the
%                           longitudinal axis used to probe the system. 
%                           Defaults to 1500.
%   rayIntersectionHeight - Scalar. Distance of the ray from the 
%                           longitudinal axis at the point the ray
%                           intersects the origin plane.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   opticalPower          - Scalar. The optical power of the system.
%   focalPoint            - 3x1 matrix. The location of the focal point.
%
% Examples:
%{
    % Determine the refractive power of the default eye
    eye = modelEyeParameters('accommodation',0);
    opticalPower = calcOpticalPower(eye);
    outline = sprintf('The refractive power of the unaccommodated model eye is %2.2f diopters.\n',opticalPower);
    fprintf(outline)
%}
%{
    % Determine the refractive power of the un-accommodated cystraline
    % lens in air
    sceneGeometry = createSceneGeometry('accommodation',0);
    opticalSystem = sceneGeometry.refraction.retinaToStop.opticalSystem;
    opticalSystem = reverseSystemDirection(opticalSystem);
    opticalPower = calcOpticalPower(opticalSystem);
    outline = sprintf('The refractive power of the unaccommodated crystaline lens is %2.2f diopters.\n',opticalPower);
    fprintf(outline)
%}

arguments
    opticalSystem
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    rayIntersectionHeight (1,1) {mustBeNumeric} = 0.5
    cameraMedium = 'air'
end

% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);
    end
end

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayOriginDistance);

% Obtain the principal point
P = calcPrincipalPoint(opticalSystem, rayOriginDistance, rayIntersectionHeight);

% Create parallel rays in the valid direction
switch systemDirection
    case 'cameraToEye'
        R1 = quadric.normalizeRay([rayOriginDistance,-1;-rayIntersectionHeight,0;0,0]);
        R2 = quadric.normalizeRay([rayOriginDistance,-1;rayIntersectionHeight,0;0,0]);
        signD = 1;
    case 'eyeToCamera'
        R1 = quadric.normalizeRay([-rayOriginDistance,1;-rayIntersectionHeight,0;0,0]);
        R2 = quadric.normalizeRay([-rayOriginDistance,1;rayIntersectionHeight,0;0,0]);
        signD = -1;
    otherwise
        error(['Not a valid system direction: ' systemDirection])
end

% Trace the rays
M1 = rayTraceQuadrics(R1, opticalSystem);
M2 = rayTraceQuadrics(R2, opticalSystem);

% Calculate the focal point from the output rays
focalPoint = quadric.distanceRays(M1,M2);

% Obtain the refractive index of the media in which the ray terminates
refractiveIndex = opticalSystem(end,end);

% Obtain the power of the optical system in diopters
opticalPower = signD * refractiveIndex / ((P(1) - focalPoint(1))/1000);

end

