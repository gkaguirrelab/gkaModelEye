function [G1, distanceError, angleError] = geodesicByReckoning( S,G0,targetDistance,targetAngle )
% Returns destination geodetic given start position, angle, and distance
%
% Syntax:
%  [G1, distanceError, angleError] = quadric.geodesicByReckoning( S,G0,targetAngle,targetDistance )
%
% Description:
%   The "direct" geodetic problem on a tri-axial ellipsoid consists of
%   finding a point (G1, specified in geodetic coordinates) that lies at a
%   specified distance and starting angle from another point (G0). This
%   routine solves the "direct" geodetic problem by an ugly, brute force
%   non-linear search over solutions of the "indirect" solution of G Panou.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   G0                    - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   targetAngle           - Scalar. Distance of the geodetic between the
%                           two points.
%   startAngle            - Scalar. The heading of the geodetic path, in
%                           degrees, relative to a line of constant omega
%                           in the ellipsoidal geodetic coordinate system
%                           on the ellipsoidal surface. 
%
% Outputs:
%   G1                    - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   distanceError         - Scalar. Difference between the provided target
%                           distance and the geodetic distance between G0
%                           and G1.
%   angleError            - Scalar. Difference in angle (in degrees)
%                           between the target and recovered start angle.
%
% Examples:
%{
    eye = modelEyeParameters('sphericalAmetropia',0);
    opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName', 'retinaToPupil' );
    S = opticalSystem(1,1:10);
    G0 = [-90;-90;0];
    G1 = [-87;-15;0];
    [distance,startAngle] = quadric.panouGeodesicDistance(S,G0,G1);
    tic
    [G1prime, distanceError, angleError] = quadric.geodesicByReckoning(S,G0,startAngle,distance);
    toc
    assert(max(abs(G1-G1prime))<1e-4);
%}

% Obtain the radii of the ellipsoid in canonical order (small to large)
radii = quadric.radii(quadric.alignAxes(S));

% For the purposes of making an initial guess as to the coordinates of G1,
% treat the ellipsoid as bi-axial, and use the MATLAB reckon function.
a = mean(radii(3)); % The semi-major axis
b = radii(1); % The semi-minor axis
x0 = reckon(G0(1),G0(2),targetDistance,targetAngle,[a sqrt(a^2-b^2)/a]);

% Define the upper and lower bounds
if x0(2) < 0
    lb = [-89.9 -180];
    ub = [89.9 0];
else
    lb = [-89.9 0];
    ub = [89.9 180];
end

% Define some nested variables
angleAtLast = [];
distanceAtLast = [];
xLast = [];

% fmincon options
options = optimoptions(@fmincon,...
    'Display','off', ...
    'Algorithm','sqp');

% Call fmincon. Note the use of anonymous functions for the objective and
% constraint, which nested below.
[G1,distanceErrorSq,~,output] = ...
    fmincon(@objfun, x0, [], [], [], [], lb, ub, @constr, options);
% Nested objective function
    function fval = objfun(x)
        % Check if computation is necessary
        if ~isequal(x,xLast)
            [distanceAtLast,angleAtLast] = quadric.panouGeodesicDistance(S,G0,x);
            xLast = x;
        end
        % Compute objective function as Euclidean distance error
        fval = (distanceAtLast-targetDistance)^2;
    end
% Nested constraint function
    function [c,ceq] = constr(x)
        if ~isequal(x,xLast)
            [distanceAtLast,angleAtLast] = quadric.panouGeodesicDistance(S,G0,x);
            xLast = x;
        end
        % c unusued
        c = [];
        % ceq:
        ceq = abs(angleAtLast-targetAngle);
    end

% Assemble the variables to return
G1 = [G1 0]'; % Need to add a zero elevation for this surface point
distanceError = sqrt(distanceErrorSq);
angleError = output.constrviolation;

end

