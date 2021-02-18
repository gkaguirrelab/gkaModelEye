function [rayPath,error] = calcSightRay(eye,G,X,stopRadius,originRayDepth,cameraMedium)
% Ray that passes through entrance pupil center to reach the retinal coord
%
% Syntax:
%  [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance)
%
% Description
%   Given an eye and a retinal coordinate, the routine identifies the ray
%   that passes through the center of the entrance pupil, and arrives at
%   the specified retinal location. If the retinal location is the fovea,
%   then this ray is the "line of sight" for the eye.
%
%   If not defined, the radius of the aperture stop is set to provide an
%   entrance pupil diameter of ~3.5 mm, which tends to produce the highest
%   degree of acuity in normal observers. The ray origin distance is
%   assumed to 500 mm unless set.
%
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
%   stopRadius            - Scalar. Radius of the aperture stop, in mm.
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
%   error                 - Scalar with the distance of ray intersection 
%                           from retinal target (in mm)
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Pick a retinal point in geodetic coordinates, a bit away from the
    % vertex
    G = [-65,-65,0];
    % Find the sight ray
    [rayPath,error] = calcSightRay(eye,G);
    % Show the optical system and sight ray
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'rayPath',rayPath,'surfaceAlpha',0.05);
    xlim([-25 10])
%}

% Code to determine the stop radius that corresponds to a pupil diameter of
% 3.5 mm. This value is used as it is found to provide peak acuity for
% normal observers.
%{
    entranceRadius = 3.5/2;
    % Prepare scene geometry and eye pose aligned with visual axis
    sceneGeometry = createSceneGeometry();
    % Obtain the pupil area in the image for the entrance radius
    % assuming no ray tracing
    sceneGeometry.refraction = [];
    pupilImage = projectModelEye([0, 0, 0, entranceRadius],sceneGeometry);
    stopArea = pupilImage(3);
    % Add the ray tracing function to the sceneGeometry
    sceneGeometry = createSceneGeometry();
    % Search across stop radii to find the value that matches the observed
    % entrance area.
    myPupilEllipse = @(radius) projectModelEye([0, 0, 0, radius],sceneGeometry);
    myArea = @(ellipseParams) ellipseParams(3);
    myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
    stopRadius = fminunc(myObj, entranceRadius);
    outline = sprintf('A 3.5mm diameter entrance pupil corresponds to a %2.2fmm stop radius\n',stopRadius);
    fprintf(outline);
%}


% Parse inputs

if nargin<2
    error('Invalid number of input arguments');
end

if nargin==2
    X = [];
    stopRadius = 1.53;
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==3
    stopRadius = 1.53;
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==4
    originRayDepth = 500;
    cameraMedium = 'air';
end

if nargin==5
    cameraMedium = 'air';
end

% Check if the X0 value is empty. If so, derive the X0
% Cartesian coordinates from the ellipsoidal geodetic coordinate.
if isempty(X)
    S = eye.retina.S;
    X = quadric.ellipsoidalGeoToCart( G, S );
end

% Make sure that X is a column vector
if all(size(X)==[1 3])
    X = X';
end

% Create a default sceneGeometry (which has a pinhole camera)
sceneGeometry = createSceneGeometry('eye',eye,'cameraMedium',cameraMedium,'surfaceSetName',{'stopToMedium','mediumToCamera','mediumToRetina'});

% Define an objective, which is the distance of the point of intersection
% of the ray in the retina from the target coordinate. The point of origin
% on the ray is set to be at the p(1), p(2) horizontal and vertical
% locations, with a depth that sets the Euclidean distance of the point
% from the cornea equal to the originRayDepth variable.
myObj = @(p) retinalDistance([sqrt(originRayDepth^2-p(1)^2-p(2)^2);p(1);p(2)],sceneGeometry,stopRadius,X);

% Define an p0 location, which is in line with the retinal target and the
% corneal apex
p0 = (originRayDepth./X(1)) .* X(2:3);

% Perform the search, finding the [h,v] location of the ray origin
p = fminsearch(myObj,p0);

% Assemble the origin point of the ray
T = [sqrt(originRayDepth^2-p(1)^2-p(2)^2);p(1);p(2)];

% Perfrom the ray trace one more time and save the output for return
[error, outputRay, rayPath] = retinalDistance(T,sceneGeometry,stopRadius, X);

% Concatenate the outputRay onto the rayPath
rayPath(:,end+1)=outputRay(:,1);


end

%% Local function
function [fVal, outputRay, rayPath] = retinalDistance(T,sceneGeometry,stopRadius, X)

% Obtain the center of the entrance pupil as seen from T
sceneGeometry.cameraPosition.translation = T([2 3 1]);
[~, ~, ~, ~, ~, eyePoints, pointLabels] = projectModelEye([0 0 0 stopRadius], sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:))';

% Trace a ray from T, towards the entrancePupilCenter, all the way to the
% retina.
R = quadric.normalizeRay([T,entrancePupilCenter-T]);
[outputRay, rayPath] = rayTraceQuadrics(R, sceneGeometry.refraction.mediumToRetina.opticalSystem);

% The error is the Euclidean distance between the retinal target coordinate
% and the intersection of the ray upon the retina.
fVal = norm(outputRay(:,1)-X);

end
