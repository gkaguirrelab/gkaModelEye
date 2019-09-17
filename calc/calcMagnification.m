function magnification = calcMagnification(opticalSystem)
% Calcuates the magnification produced by an optical system
%
% Syntax:
%  magnification = calcDiopters(opticalSystem)
%
% Description
%   Calculates the angular magnification produced by an optical system.
%   Values greater than 1 indicate magnification, values less than one
%   reflect minification.
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
%   magnification         - Scalar. The magnification produced by the
%                           system.
%
% Examples:
%{
    % Create a lens and report the magnification it produces
    lensDiopters = 5;
    opticalSystem=addSpectacleLens([],lensDiopters);
    magnification = calcMagnification(opticalSystem)
%}


% Obtain the systemDirection
systemDirection = calcSystemDirection(opticalSystem);

% Create a ray in the valid direction
switch systemDirection
    case 'cameraToEye'
        angleInitial = 179;
        R = quadric.normalizeRay(quadric.anglesToRay([50;0;0], angleInitial, 0 ));
    case 'eyeToCamera'
        angleInitial = 1;
        R = quadric.normalizeRay(quadric.anglesToRay([-3.9;0;0], angleInitial, 0 ));
    otherwise
        error('Not a valid system direction')
end

% Trace a ray through the system
outputRay = rayTraceQuadrics(R, opticalSystem);

% Obtain the angle of the output ray w.r.t. the optical axis
angleFinal = quadric.rayToAngles(outputRay);
magnification = angleInitial / angleFinal;
    
end

