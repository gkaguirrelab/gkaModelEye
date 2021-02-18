function [rayPath,nodalPoints,errors] = calcNodalRay(eye,G,X,originRayDepth,cameraMedium)
% Returns the path of the nodal ray that intersects the retinal coordinate
%
% Syntax:
%  [rayPath,nodes,errors] = calcNodalRay(eye,G,X,originRayDepth,cameraMedium)
%
% Description
%   Given an eye structure and a coordinate, the routine returns a matrix
%   that contains the path of a ray that arrives at this coord and has an
%   angle of incidence at cornea (w.r.t the optical axis) equal to the
%   angle with which it intersects the retina. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The routine returns the nodalPoints,which are found by extending the
%   initial and exit segments of the ray to the optical axis. The
%   nodal ray that arrives at the fovea is termed the "visual axis" of the
%   eye.
%
%   For a model eye with decentered and/or astigmatic elements, there may
%   not exist a true nodal ray, but only approximations to this concept.
%   There may not be a ray that both has incident and emergent rays, and
%   intersects the retinal target. Further, the incident and emergent rays,
%   when extended, will pass by but not intersect the optical axis. These
%   various error are supplied in the errors output variable. Further, the
%   particular nodal points will vary depending upon the particular retinal
%   target.
%
%   The routine can accept points on the ellipsoidal surface of the retina
%   specified in either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   G                     - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X                     - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   originRayDepth        - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
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
%                             - distance of ray intersection from retinal
%                               target (in mm)
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
    % Pick a retinal point in geodetic coordinates, a bit away from the
    % vertex
    G = [-65,-65,0];
    % Find the nodal ray
    [rayPath,nodalPoints,errors] = calcNodalRay(eye,G);
    % Show the optical system, nodal ray, and nodal points
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'rayPath',rayPath,'surfaceAlpha',0.05);
    hold on
    xlim([-25 10])
    plot3(nodalPoints(1,:),nodalPoints(2,:),nodalPoints(3,:),'*b')
%}


% Parse inputs
if nargin<2
    error('Invalid number of input arguments');
end

if nargin==2
    X = [];
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==3
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
end

% Check if the X value is empty. If so, derive the X
% Cartesian coordinates from the ellipsoidal geodetic coordinate.
if isempty(X)
    S = eye.retina.S;
    X = quadric.ellipsoidalGeoToCart( G, S );
end

% Make sure that X is a column vector
if all(size(X)==[1 3])
    X = X';
end

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);


% Initialize an anonymous function for the full objective. This combines
% both the objective and the non-linear constraint.
fullObj = @(x) fullObjective(x,opticalSystem,originRayDepth,X);

% Define some search options for fminsearch
options = optimset('fmincon');
options.FunValCheck = 'off'; % Tolerate nans from the objective
options.Display = 'off'; % Silencio

% x0, and bounds. x0 is placed so that it is in line with the un-rotated
% corneal apex and the retinal target location,, with a slight jitter to
% avoid local minima in a rotationally symmetric, aligned optical system.
% The params are:
%   - horizontal position of the ray origin (in mm)
%   - vertical position of the ray origin (in mm)
%   - p1 angle of the initial ray
%   - p2 angle of the initial ray
T = (originRayDepth./X(1)) .* X;
[p(1),p(2)] = quadric.rayToAngles(quadric.normalizeRay([T,-T]));
x0 = [T(2),T(3),p(1)+0.01,p(2)-0.01];
lb = [-originRayDepth*10,-originRayDepth*10,-360,-360];
ub = [originRayDepth*10,originRayDepth*10,360,360];

% Initialize nested variables
xLast = [];
cLast = [];
ceqLast = [];
fValLast = [];
ceqFirstFlag = true;

% Perform the search, using the nested function
x = fmincon(@objFun, x0,[],[],[],[],lb,ub,@nonlcon,options);
    function fVal = objFun(x)
        if ~isequal(x,xLast)
            [fValLast,cLast,ceqLast] = fullObj(x);
            xLast = x;
        end
        fVal = fValLast;
        % This is some business to prevent fmincon from considering as
        % acceptable a ceq that is low relative to the x0 ceq
        if ceqFirstFlag
            ceqFirstFlag = false;
        else
            if ceqLast == 1e3
                ceqLast = 1e6;
            end
        end
    end
    function [c, ceq] = nonlcon(x)
        if ~isequal(x,xLast)
            [fValLast,cLast,ceqLast] = fullObj(x);
            xLast = x;
        end
        c = cLast;
        % This is some business to prevent fmincon from considering as
        % acceptable a ceq that is low relative to the x0 ceq
        ceq = ceqLast;
        if ceqFirstFlag
            ceqFirstFlag = false;
        else
            if ceq == 1e3
                ceq = 1e6;
            end
        end
    end

% Evaluate the objective function once more, using the final values
[retinalDistanceError, ~, angleError, outputRay,rayPath] = ...
    fullObjective(x,opticalSystem,originRayDepth,X);

% Find the nodal points
opticalAxis = [0,1;0,0;0,0];
inputRay = quadric.normalizeRay([rayPath(:,1),rayPath(:,2)-rayPath(:,1)]);
[~,iNodeError,incidentNode] = quadric.distanceRays(inputRay,opticalAxis);
[~,eNodeError,emergentNode] = quadric.distanceRays(outputRay,opticalAxis);

% Assemble the errors and nodal points for return
nodalPoints = [incidentNode,emergentNode];
errors = [retinalDistanceError,angleError,iNodeError,eNodeError];

% Concatenate the outputRay onto the rayPath
rayPath(:,end+1)=outputRay(:,1);

end


%% Local function

function [fVal, c, ceq, outputRay,rayPath] = fullObjective(x,opticalSystem,originRayDepth,X)

% Define a point in space
T = [sqrt(originRayDepth^2-x(1)^2-x(2)^2);x(1);x(2)];

% Trace from T at angle p.
R = quadric.anglesToRay(T,x(3), x(4));
[outputRay,rayPath] = rayTraceQuadrics(R, opticalSystem);

% c is unused
c = [];

% Find the difference in angles between the ray leaving T, and the ray
% arriving at the retina
ceq = quadric.angleRays( R, outputRay );

% Find the Euclidean distance between the retinal target coordinate and the
% intersection of the ray upon the retina.
fVal = norm(outputRay(:,1)-X);

end