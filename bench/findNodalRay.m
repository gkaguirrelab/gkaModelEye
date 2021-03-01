function [rayPath,angleError] = findNodalRay(rayOrigin,opticalSystem,incidentNodeX0)
% Returns the path of the nodal ray from the starting coordinate
%
% Syntax:
%  [rayPath,errors] = findNodalRay(rayOrigin,opticalSystem,incidentNodeX0)
%
% Description
%   Given an opticalSystem and a coordinate in eye coordinate space, the
%   routine returns a matrix that contains the path of a ray that departs
%   from this coord and has an angle of incidence at the first optical
%   surface equal to the angle with which it leaves the last surface. This
%   is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
% Inputs:
%   rayOrigin             - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in the eye coordinate space with the
%                           dimensions [p1, p2, p3] from which the nodal
%                           ray will originate.
%   opticalSystem         - Struct. See assembleOpticalSystem.m
%   incidentNodeX0        - An optional 1x3 vector that gives the location
%                           in eye  space that is an initial guess for the
%                           location of the incident node of the optical
%                           system. If not supplied, a value that is
%                           typical for the human eye is used.
%
% Outputs:
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   angleError            - Scalar. The departure from parallel of the 
%                           incident and emergent rays (deg)
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Obtain the optical system
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    % Define a point in eye-world coordinate space
    X = [150, 20, 10];
    % Find the nodal ray
    [rayPath,angleError] = findNodalRay(X,opticalSystem);
%}

if nargin==2
    incidentNodeX0 = [-7 0 0];
end

% Place the vectors in column orientation
rayOrigin = rayOrigin';
incidentNodeX0 = incidentNodeX0';

% Set an guess for the angles of the initial ray
p0 = nan(1,2);
[p0(1),p0(2)] = quadric.rayToAngles(quadric.normalizeRay([rayOrigin,incidentNodeX0-rayOrigin]));

% Initialize an anonymous function for the objective
myObj = @(p) objective(p,rayOrigin,opticalSystem);

% Search
p = fminsearch(myObj,p0);

% Evaluate the objective function once more, using the found values
[angleError,rayPath] = objective(p,rayOrigin,opticalSystem);

end


%% Local function

function [fVal,rayPath] = objective(p,rayOrigin,opticalSystem)

% Trace from rayOrigin at p angles.
R = quadric.anglesToRay(rayOrigin,p(1),p(2));
[outputRay,rayPath] = rayTraceQuadrics(R, opticalSystem);

% Find the absolute difference in angles between the ray leaving T, and the
% ray arriving at the retina
fVal = abs(quadric.angleRays( R, outputRay ));

end