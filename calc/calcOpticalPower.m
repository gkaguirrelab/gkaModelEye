function [diopters, focalPoint] = calcOpticalPower(opticalSystem, rayStartDepth, rayHeight)
% Calcuate the refractive power of an opticalSystem
%
% Syntax:
%  [diopters, focalPoint] = calcOpticalPower(opticalSystem)
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
%   Note that for optical systems with spherical aberration, the calculated
%   optical power will vary depending upon the path of the ray. The
%   rayHeight parameter controls the height at which the ray strikes the
%   first surface.
%
% Inputs:
%   opticalSystem         - An mx19 matrix, where m is set by the key value
%                           opticalSystemNumRows. Each row contains the
%                           values:
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
%   rayStartDepth         - Scalar. Point of origin of the ray along the
%                           optical axis used to probe the system. Defaults
%                           to 100.
%   rayHeight             - Scalar. Distance of the ray origin from the
%                           optical axis.
%
% Outputs:
%   diopters              - Scalar. The optical power of the system.
%   focalPoint            - 3x1 matrix. The location of the focal point.
%
% Examples:
%{
    % Determine the refractive power of the unaccommodated eye
    sceneGeometry = createSceneGeometry();
    diopters = calcOpticalPower(sceneGeometry.refraction.cameraToRetina.opticalSystem);
    outline = sprintf('The refractive power of the unaccommodated model eye is %2.2f diopters.\n',diopters);
    fprintf(outline)
%}
%{
    % Determine the refractive power of the un-accommodated cystraline
    % lens in air
    sceneGeometry = createSceneGeometry();
    opticalSystem = sceneGeometry.refraction.retinaToStop.opticalSystem;
    opticalSystem = reverseSystemDirection(opticalSystem);
    diopters = calcOpticalPower(opticalSystem);
    outline = sprintf('The refractive power of the unaccommodated crystaline lens is %2.2f diopters.\n',diopters);
    fprintf(outline)
%}

% Handle nargin
if nargin==1
    rayStartDepth = 100;
    rayHeight = 1;
end

% Handle nargin
if nargin==2
    rayHeight = 1;
end


% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayStartDepth);

% Obtain the principal point
P = calcPrincipalPoint(opticalSystem, rayStartDepth, rayHeight);

% Create parallel rays in the valid direction
switch systemDirection
    case 'cameraToEye'
        R1 = quadric.normalizeRay([rayStartDepth,-1;-rayHeight,0;0,0]);
        R2 = quadric.normalizeRay([rayStartDepth,-1;rayHeight,0;0,0]);
        signD = 1;
    case 'eyeToCamera'
        R1 = quadric.normalizeRay([-rayStartDepth,1;-rayHeight,0;0,0]);
        R2 = quadric.normalizeRay([-rayStartDepth,1;rayHeight,0;0,0]);
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
diopters = signD * refractiveIndex / ((P(1) - focalPoint(1))/1000);

end

