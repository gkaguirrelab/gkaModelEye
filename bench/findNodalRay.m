function [rayPath,nodalPoints,errors] = findNodalRay(eyeCoordOrigin,opticalSystem,incidentNodeX0)
% Returns the path of the nodal ray from the starting coordinate
%
% Syntax:
%  [rayPath,nodalPoints,errors] = findNodalRay(eyeCoordOrigin,opticalSystem,incidentNodeX0)
%
% Description
%   Given an opticalSystem and a coordinate in eye coordinate space, the
%   routine returns a matrix that contains the path of a ray that departs
%   from this coord and has an angle of incidence at the first optical
%   surface (w.r.t the optical axis) equal to the angle with which it
%   leaves the last surface. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The routine returns the nodalPoints,which are found by extending the
%   initial and exit segments of the ray to the optical axis.
%
% Inputs:
%   eyeCoordOrigin        - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in the eye coordinate space with the
%                           dimensions p1, p2, p3.
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
%   nodalPoints           - 3x2 matrix that provides the approximation to
%                           incident and emergent nodal points in the eye
%                           coordinate space. This is the point on each ray
%                           that is closest to the optical axis.
%   errors                - 1x4 matrix with the follow error values:
%                             - departure from parallel of the incident and
%                               emergent rays (deg)
%                             - distance of the incident nodal point from
%                               the incident ray
%                             - distance of the emergent nodal point from
%                               the emergent ray
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
    [rayPath,nodalPoints,errors] = findNodalRay(X,opticalSystem);
    % Show the optical system, nodal ray, and nodal points
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'rayPath',rayPath,'surfaceAlpha',0.05);
    hold on
    xlim([-25 10])
    plot3(nodalPoints(1,:),nodalPoints(2,:),nodalPoints(3,:),'*b')
%}

if nargin==2
    incidentNodeX0 = [-7 0 0];
end

% Place the vectors in column orientation
eyeCoordOrigin = eyeCoordOrigin';
incidentNodeX0 = incidentNodeX0';

% Set an guess for the angles of the initial ray
p0 = nan(1,2);
[p0(1),p0(2)] = quadric.rayToAngles(quadric.normalizeRay([eyeCoordOrigin,incidentNodeX0-eyeCoordOrigin]));

% Initialize an anonymous function for the objective
myObj = @(p) objective(p,eyeCoordOrigin,opticalSystem);

% Search
p = fminsearch(myObj,p0);

% Evaluate the objective function once more, using the found values
[angleError,outputRay,rayPath] = objective(p,eyeCoordOrigin,opticalSystem);

% Find the nodal points
opticalAxis = [0,1;0,0;0,0];
inputRay = quadric.normalizeRay([rayPath(:,1),rayPath(:,2)-rayPath(:,1)]);
[~,iNodeError,incidentNode] = quadric.distanceRays(inputRay,opticalAxis);
[~,eNodeError,emergentNode] = quadric.distanceRays(outputRay,opticalAxis);

% Assemble the errors and nodal points for return
nodalPoints = [incidentNode,emergentNode];
errors = [angleError,iNodeError,eNodeError];

% Concatenate the outputRay onto the rayPath
rayPathFull = nan(3,size(rayPath,2)+1);
rayPathFull(:,1:end-1)=rayPath;
rayPathFull(:,end)=outputRay(:,1);
rayPath = rayPathFull;

end


%% Local function

function [fVal, outputRay,rayPath] = objective(p,eyeCoordOrigin,opticalSystem)


% Trace from eyeCoordOrigin at p angles.
R = quadric.anglesToRay(eyeCoordOrigin,p(1),p(2));
[outputRay,rayPath] = rayTraceQuadrics(R, opticalSystem);

% Find the absolute difference in angles between the ray leaving T, and the
% ray arriving at the retina
fVal = abs(quadric.angleRays( R, outputRay ));


end