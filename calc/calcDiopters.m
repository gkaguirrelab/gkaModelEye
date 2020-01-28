function [diopters, focalPoint] = calcDiopters(opticalSystem, forceEyeToCamera, rayStartDepth)
% Calcuate the cameraToEye-direction refractive power of an opticalSystem
%
% Syntax:
%  [diopters, focalPoint] = calcDiopters(opticalSystem)
%
% Description
%   Calculates the refractive power of an optical system in units of
%	diopters. A negative value specifies a lens that would be worn by
%	someone with myopia to correct their vision.
%
%   Some optical systems end in a medium with a refractive index other than
%   one. In this case the diopter strength is given by:
%
%       diopters = refractiveIndex / effectiveFocalLength
%
%   where the effective focal length is the distance between the principal
%   point of the optical system and the focal point.
%
%   This routine returns optical power in the cameraToEye direction for the
%   passed optical system, unless the "forceEyeToCamera" argument is set to
%   true.
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
%   forceEyeToCamera      - Logical. Optional. Defaults to false if not
%                           set.
%   rayStartDepth         - Scalar. Point of origin of the ray along the
%                           optical axis used to probe the system. Defaults
%                           to 100.
%
% Outputs:
%   diopters              - Scalar. The optical power of the system.
%                           Calculated in the cameraToEye direction unless
%                           forceEyeToCamera is set to true.
%   focalPoint            - 3x1 matrix. The location of the focal point.
%
% Examples:
%{
    % Determine the refractive power of the model eye
    eye=modelEyeParameters('navarroD',0);
    opticalSystem=assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[]);
    [diopters, focalPoint] = calcDiopters(opticalSystem)
%}
%{
    % Determine the refractive power of the lens 
    eye=modelEyeParameters('accommodationDiopters',0);
    opticalSystem=assembleOpticalSystem(eye,'surfaceSetName','stopToRetina','opticalSystemNumRows',[]);
    [diopters, focalPoint] = calcDiopters(opticalSystem)
%}

% Handle nargin
if nargin==1
    forceEyeToCamera = false;
    rayStartDepth = [100, -100];
end

% Handle nargin
if nargin==2
    rayStartDepth = [100, -100];
end

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayStartDepth);

% Unless the forceEyeToCamera flag is set, ensure that the optical system
% is in the cameraToEye state.
if strcmp(systemDirection,'eyeToCamera') && ~forceEyeToCamera
    opticalSystem = reverseSystemDirection(opticalSystem);
end
if strcmp(systemDirection,'cameraToEye') && forceEyeToCamera
    opticalSystem = reverseSystemDirection(opticalSystem);
end

% Obtain the system direction again after that potential reversing
systemDirection = calcSystemDirection(opticalSystem, rayStartDepth);

% Obtain the principal point
P = calcPrincipalPoint(opticalSystem, rayStartDepth);

% Create parallel rays in the valid direction
switch systemDirection
    case 'cameraToEye'
        R1 = quadric.normalizeRay([rayStartDepth(1),-1;-1,0;0,0]);
        R2 = quadric.normalizeRay([rayStartDepth(1),-1;1,0;0,0]);
        signD = 1;
    case 'eyeToCamera'
        R1 = quadric.normalizeRay([rayStartDepth(2),1;-1,0;0,0]);
        R2 = quadric.normalizeRay([rayStartDepth(2),1;1,0;0,0]);
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

