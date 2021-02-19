function [rayPath,error] = calcSightRayToRetina(eye,rayDestination,rayOriginDistance,stopRadius,cameraMedium)
% Ray that passes through entrance pupil center to reach the retinal coord
%
% Syntax:
%  [outputRay,rayPath,fixationEyePose,foveaDistanceError] = calcSightRayToRetina(sceneGeometry,stopRadius,fixTargetDistance)
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
%   rayDestination        - 1x3 or 3x1 vector that provides the coordinates
%                           of a location on the retina. The location can
%                           be given in Cartesian or geodetic coordinates.
%                           For the latter, the values are beta, omega, and
%                           elevation in units of degrees. Beta is defined
%                           over the range -90:90, and omega over the range
%                           -180:180. Elevation has an obligatory value of
%                           zero as this solution is only defined on the
%                           surface. The routine will attempt to determine
%                           which type of coordinate set has been provided.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the corneal apex. Assumed to be
%                           500 mm if not defined.
%   stopRadius            - Scalar. Radius of the aperture stop, in mm.
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
    % Obtain the coordinates of the fovea
    rayDestination = eye.landmarks.fovea.coords;
    % Find the sight ray to the fovea (i.e., the line of sight axis)
    [rayPath,error] = calcSightRayToRetina(eye,rayDestination);
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
    error('calcSightRayToRetina:invalidArguments','Too few input arguments')
end

if nargin==2
    rayOriginDistance = 500;
    stopRadius = 1.53;
    cameraMedium = 'air';
end

if nargin==3
    stopRadius = 1.53;
    cameraMedium = 'air';
end

if nargin==4
    cameraMedium = 'air';
end

% Make rayDestination a column vector
if all(size(rayDestination)==[1 3])
    rayDestination = rayDestination';
end

% Time to interpret the rayDestination variable. Hard-code a tolerance for
% the distance of the rayDestination coordinate from the retinal surface
surfaceTol = 1e-6;

% Get the retinal surface and quadric function
S = eye.retina.S;
Sfunc = quadric.vecToFunc(S);

% Interpret the value as Cartesian, if it is not a valid retinal location,
% the quadric function will evaluated with a non-zero value
if Sfunc(rayDestination(1),rayDestination(2),rayDestination(3)) > surfaceTol
    
    % If the last value of the coordinate is not zero, then this can't be a
    % geodetic coordinate either
    if rayDestination(3) > surfaceTol
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Check that the candidate beta and omega values are in range
    if abs(rayDestination(1))>90 || abs(rayDestination(2))>180
        error('calcSightRayToRetina:invalidCoordinate','Supply a Cartesian or geodetic coordinate that is on the retinal surface.')
    end
    
    % Looks like a valid geodetic coordinate. Convert to Cartesian.
    rayDestination = quadric.ellipsoidalGeoToCart( rayDestination, S );
        
end

% Create a default sceneGeometry (which has a pinhole camera)
sceneGeometry = createSceneGeometry('eye',eye,...
    'cameraMedium',cameraMedium,...
    'surfaceSetName',{'stopToMedium','mediumToCamera','mediumToRetina'});

% Define an objective, which is the distance of the point of intersection
% of the ray in the retina from the target coordinate. The point of origin
% on the ray is set to be at the p(1), p(2) horizontal and vertical
% locations, with a depth that sets the Euclidean distance of the point
% from the cornea equal to the rayOriginDistance variable.
myCoord = @(p) [sqrt(rayOriginDistance^2-p(1)^2-p(2)^2);p(1);p(2)];
myObj = @(p) objective(myCoord(p),sceneGeometry,stopRadius,rayDestination);

% Define an p0 location, which is in line with the retinal target and the
% corneal apex
p0 = (rayOriginDistance./rayDestination(1)) .* rayDestination(2:3);

% Perform the search, finding the [h,v] location of the ray origin
p = fminsearch(myObj,p0);

% Assemble the origin point of the ray
rayOrigin = myCoord(p);

% Perfrom the ray trace one more time and save the output for return
[error, outputRay, rayPath] = objective(rayOrigin,sceneGeometry,stopRadius,rayDestination);

% Concatenate the outputRay onto the rayPath
rayPath(:,end+1)=outputRay(:,1);


end


%% Local function
function [fVal, outputRay, rayPath] = objective(rayOrigin,sceneGeometry,stopRadius, rayDestination)

% Obtain the center of the entrance pupil as seen from T
sceneGeometry.cameraPosition.translation = rayOrigin([2 3 1]);
[~, ~, ~, ~, ~, eyePoints, pointLabels] = projectModelEye([0 0 0 stopRadius], sceneGeometry, 'nStopPerimPoints', 16);
entrancePupilCenter = mean(eyePoints(strcmp(pointLabels,'pupilPerimeter'),:))';

% Trace a ray from T, towards the entrancePupilCenter, all the way to the
% retina.
R = quadric.normalizeRay([rayOrigin,entrancePupilCenter-rayOrigin]);
[outputRay, rayPath] = rayTraceQuadrics(R, sceneGeometry.refraction.mediumToRetina.opticalSystem);

% The error is the Euclidean distance between the retinal target coordinate
% and the intersection of the ray upon the retina.
fVal = norm(outputRay(:,1)-rayDestination);

end
