function [diopters, focalPoint] = calcDiopters(opticalSystem)
% Calcuates the optical power of an opticalSystem
%
% Syntax:
%  [diopters, focalPoint] = calcDiopters(opticalSystem)
%
% Description
%   Calculates the power refractive power of an optical system in units of
%	diopters. A negative value specifies a lens that would be worn by
%	someone with myopia to correct their vision.
%
%   The implementation of optical systems and ray tracing in this code
%   results in only one of these ray tracing directions being available for
%   a given opticalSystem variable. We try both directions here and report
%   back the valid solution.
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
%
% Outputs:
%   diopters              - Scalar. The optical power of the system.
%   focalPoint            - 3x1 matrix. The location of the focal point.
%
% Examples:
%{
    % Display a lens and it's principal point
    opticalSystem = addSpectacleLens([],-5);
    % Plot this
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true);
    % Find the principal points and plot these
    [P, systemDirectionOut] = calcPrincipalPoints(opticalSystem);
    plot3(P(1),P(2),P(3),'*r')
    assert(strcmp(systemDirection,systemDirectionOut));
%}

% Obtain the principal point and systemDirection
P = calcPrincipalPoint(opticalSystem);
systemDirection = calcSystemDirection(opticalSystem);

% Create parallel rays in the valid direction
switch systemDirection
    case 'cameraToEye'
        R1 = quadric.normalizeRay([100,-1;-3,0;0,0]);
        R2 = quadric.normalizeRay([100,-1;3,0;0,0]);
        signD = 1;
    case 'eyeToCamera'
        R1 = quadric.normalizeRay([-100,1;-3,0;0,0]);
        R2 = quadric.normalizeRay([-100,1;3,0;0,0]);
        signD = -1;
    otherwise
        error('Not a valid system direction')
end

% Trace the rays
M1 = rayTraceQuadrics(R1, opticalSystem);
M2 = rayTraceQuadrics(R2, opticalSystem);

% Calculate the focal point from the output rays
focalPoint = quadric.distanceRays(M1,M2);

% Obtain the power of the optical system in diopters
diopters = signD * 1000 / (P(1) - focalPoint(1));


end

