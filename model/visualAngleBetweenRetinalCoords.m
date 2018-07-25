function visualAngles = visualAngleBetweenRetinalCoords(sceneGeometry,G0,G1,X0,X1)
% The visual angles between two retinal points
%
% Syntax:
%  visualAngle = visualAngleBetweenRetinalCoords(sceneGeometry,G0,G1,X0,X1)
%
% Description
%   Given a sceneGeometry and two coordinates on the retinal surface, the
%   routine returns a vector that contains the visual angle (in degrees) 
%   between the two points, projected on the p1p2 and p1p3 planes (i.e.,
%   horizontal and vertical visual angle).
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   G0, G1                - 3x1 vectors that provide the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%
% Outputs:
%   visualAngles          - 1x2 vector with the visual angle, in degrees 
%                           between the two points within the p1p2 and p1p3
%                           planes.
%
% Examples:
%{
%}

% If only three input values were passed, derive the X0/X1 Cartesian
% coordinates from the ellipsoidal geodetic coordinates.
if nargin==3
    S = sceneGeometry.refraction.retinaToPupil.opticalSystem(1,1:10);
    X0 = quadric.ellipsoidalGeoToCart( G0, S )';
    X1 = quadric.ellipsoidalGeoToCart( G1, S )';
end

% Assemble arguments for the virtual image function to trace from the
% retina to the center of the pupil aperture. The "camera" position is set
% as the pupil center, with the dimensions rearranged for world
% coordinates.
args = {sceneGeometry.eye.pupil.center([2 3 1])', ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.refraction.retinaToPupil.opticalSystem};

% Set eyePose to all zeros (i.e., the eye is not rotated).
eyePose = [0 0 0 0];

% Ray trace to the center of the pupil
R0 = sceneGeometry.refraction.handle(X0, eyePose, args{:});
R1 = sceneGeometry.refraction.handle(X1, eyePose, args{:});

% Normalize the rays, and reverse the direction so that the ray is headed
% from the pupil center out towards the cornea
R0 = quadric.normalizeRay(R0');
R0(:,2)=-R0(:,2);
R1 = quadric.normalizeRay(R1');
R1(:,2)=-R1(:,2);

% Ray trace from the pupil center through the cornea
R0 = rayTraceQuadrics(R0, sceneGeometry.refraction.pupilToCamera.opticalSystem);
R1 = rayTraceQuadrics(R1, sceneGeometry.refraction.pupilToCamera.opticalSystem);

% Calculate and return the signed angles between the two rays
[~, angle_p1p2, angle_p1p3] = quadric.angleRays( R0, R1 );
visualAngles = [angle_p1p2, angle_p1p3];

end