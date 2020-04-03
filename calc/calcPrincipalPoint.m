function P = calcPrincipalPoint(opticalSystem, rayStartDepth)
% Returns the principal point for an opticalSystem
%
% Syntax:
%  P = calcPrincipalPoint(opticalSystem, rayStartDepth)
%
% Description
%   For a ray that passes through an optical system, the incoming and
%   outgoing rays may be extrapolated until they intersect. The point of
%   intersection is termed the "equivalent refracting locus". A set of
%   parallel, incoming rays define a corresponding set of loci, which in
%   turn define the "principal surface". The point at which this surface
%   intersects the optical axis of the system is termed a principal point
%   (P).
%
%   An optical system has two principal points, which correspond to rays
%   originating from the right or left. The implementation of optical
%   systems and ray tracing in this code results in only one of these
%   ray tracing directions being available for a given opticalSystem
%   variable; only the valid solution is returned.
%
%   In the paraxial approximation, nodal points and principal points are
%   the same.
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
%   P                     - 3x1 matrix. The location of the principal
%                           point.
%
% Examples:
%{
    % Display a lens and it's principal points
    systemDirection = 'cameraToEye';
    opticalSystem = addSpectacleLens([],-10,'systemDirection',systemDirection);
    % Plot this
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true);
    % Find the principal point and plot this
    P = calcPrincipalPoint(opticalSystem);
    plot3(P(1),P(2),P(3),'*r')
%}

% Handle nargin
if nargin==1
    rayStartDepth = [100, -100];
end

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayStartDepth);

% Create a paraxial ray in the valid direction
switch systemDirection
    case 'cameraToEye'
        R = quadric.normalizeRay(quadric.anglesToRay([rayStartDepth(1);0;0],180-0.01,0));
    case 'eyeToCamera'
        R = quadric.normalizeRay(quadric.anglesToRay([rayStartDepth(2);0;0],0.01,0));
    otherwise
        error(['Not a valid system direction: ' systemDirection])
end

% Trace the ray and find the intersection of the incoming and outgoing ray 
M = rayTraceQuadrics(R, opticalSystem);
P = quadric.distanceRays(R,M);

% The first element of P is the position along the optical axis. The other
% values are set to zero as we are working with a paraxial approximation.
P(2:3)=0;

end

