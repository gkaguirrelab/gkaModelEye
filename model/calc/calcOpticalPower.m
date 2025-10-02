function [opticalPower, focalPoint] = calcOpticalPower(opticalSystem, rayOriginDistance, rayIntersectionHeight, cameraMedium)
% Calcuate the refractive power of an opticalSystem
%
% Syntax:
%  [opticalPower, focalPoint] = calcOpticalPower(opticalSystem, rayOriginDistance, rayIntersectionHeight, cameraMedium)
%
% Description
%   Calculates the refractive power of an optical system in units of
%   diopters. A negative value specifies the power of a lens that would be
%   worn by someone with myopia to correct their vision.
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
%   opticalPower          - 3x1 matrix corresponding to sphere, cylinder, 
%                           and axis values for the system, following the
%                           "Plus Cylinder" convention, and the convention
%                           that 0/180 degrees for the cylinder corresponds
%                           to the horizontal meridian.
%   focalPoint            - 3x1 matrix. The location of the focal point.
%
% Examples:
%{
    % Determine the refractive power of the default eye
    eye = modelEyeParameters('accommodation',0);
    opticalPower = calcOpticalPower(eye);
    outline = sprintf('The sphere, cylinder of the unaccommodated model eye is %2.2f, %2.2f diopters, and the axis is %2.1f degrees.\n',opticalPower);
    fprintf(outline)
%}
%{
    % Determine the refractive power of the un-accommodated crystaline
    % lens in air
    sceneGeometry = createSceneGeometry('accommodation',0);
    opticalSystem = sceneGeometry.refraction.retinaToStop.opticalSystem;
    opticalSystem = reverseSystemDirection(opticalSystem);
    opticalPower = calcOpticalPower(opticalSystem);
    outline = sprintf('The refractive power of the unaccommodated crystaline lens is %2.2f diopters.\n',opticalPower(1));
    fprintf(outline)
%}

arguments
    opticalSystem {mustBeOpticalSystemCapable}
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    rayIntersectionHeight (1,1) {mustBeNumeric} = 0.5
    cameraMedium = 'air'
end

% Create the optical system
opticalSystem = parseOpticalSystemArgument(opticalSystem,'mediumToRetina',cameraMedium);

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayOriginDistance);

% Obtain the principal point
P = calcPrincipalPoint(opticalSystem, rayOriginDistance, rayIntersectionHeight);

% Create parallel rays in the valid direction across a variety of
% orientations. This allows us to measure cylinder and axis
for pp = 0:179
    switch systemDirection
        case 'cameraToEye'
            R1 = quadric.normalizeRay([rayOriginDistance,-1;sind(pp)*rayIntersectionHeight,0;cosd(pp)*rayIntersectionHeight,0]);
            R2 = quadric.normalizeRay([rayOriginDistance,-1;sind(pp+180)*rayIntersectionHeight,0;cosd(pp+180)*rayIntersectionHeight,0]);
            signD = 1;
        case 'eyeToCamera'
            R1 = quadric.normalizeRay([-rayOriginDistance,1;sind(pp)*rayIntersectionHeight,0;cosd(pp)*rayIntersectionHeight,0]);
            R2 = quadric.normalizeRay([-rayOriginDistance,1;sind(pp+180)*rayIntersectionHeight,0;cosd(pp+180)*rayIntersectionHeight,0]);
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
    opticalPowerByAngle(pp+1) = signD * refractiveIndex / ((P(1) - focalPoint(1))/1000);
end

% Report optical power, cylinder, and axis using the "Plus Cylinder"
% convention.
[sphereDiopters,axisDegrees] = min(opticalPowerByAngle);
axisDegrees = axisDegrees-1;
cylinderDiopters = max(opticalPowerByAngle)-min(opticalPowerByAngle);

opticalPower = [sphereDiopters,cylinderDiopters,axisDegrees];

end

