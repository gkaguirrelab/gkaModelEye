function [entranceWindowCenter,objectCoord,entranceWindowPerimeter,radiusEntranceWindow] = calcEntranceWindow(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,nStopPerimPoints,paraxialThresh,cameraMedium)
% Returns the coordinates of the entrance window of an eye
%
% Syntax:
%  [entranceWindowCenter,objectCoord,entranceWindowPerimeter] = calcEntranceWindow(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,nStopPerimPoints,paraxialThresh,cameraMedium)
%
% Description
%   The entrance window is the image of the field stop of an optical system
%   as viewed from an off-axis location. The field stop is the limiting
%   aperture for rays traveling from the object to the image. For the eye,
%   the aperture of the iris is the field stop. This routine returns the
%   center (and perimeter) of the entrance window of the eye as viewed from
%   a specified position, and for an aperture stop of a specified radius.
%
%   I position the entrance window so that it is centered on the
%   longitudinal axis in the horizontal and vertical directions. When the
%   object position is on the longitudinal axis, the entrance window is the
%   proper entrance pupil.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   stopRadius            - Scalar. The size in mm of the aperture stop.
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in angleReferenceCoord.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the distanceReferenceCoord.
%   angleReferenceCoord   - 3x1 vector that provides the coordinate from
%                           which the ray origin angles and distance are
%                           to be calculated. By default, this is [0;0;0],
%                           which is the origin coordinate on the
%                           longitudinal axis.
%   distanceReferenceCoord - 3x1 vector that provides the coordinate from
%                           which the rayOriginDistance is calculated. The
%                           The principal point is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   nStopPerimPoints      - Scalar. Then number of points on the stop
%                           perimeter to be evaluated.
%   paraxialThresh        - Scalar. The fieldRay will be forced to be
%                           separated from the longitudinal axis by at
%                           least this angle.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   entranceWindowCenter  - 3x1 vector that specifies the center of the
%                           entrance window when viewed from the specified
%                           field position.
%   objectCoord           - 3x1 vector that specifies the location of the
%                           object point for which the entrance window was
%                           calculated
%   entranceWindowPerimeter - 3xnStopPerimPoints. The set of points on the
%                           entrance window perimeter
%   radiusEntranceWindow  - Scalar. The radius of the entrance window.
%
% Examples:
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Find the center of the entrance window
    windowCoord = calcEntranceWindow(eye);
%}
%{
    % Plot the location of the entrance window in an optical system
    eye = modelEyeParameters();
    [~,~,perimeter] = calcEntranceWindow(eye,[30 0]);
    plotOpticalSystem(eye);
    plot3(perimeter(1,:),perimeter(2,:),perimeter(3,:),'*k');
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0,0]
    stopRadius (1,1) {mustBeNumeric} = 2
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    nStopPerimPoints (1,1) {mustBeNumeric} = 6
    paraxialThresh (1,1) {mustBeNumeric} = 1e-3
    cameraMedium = 'air'
end


% Check if we have a compiled version of findPupilRay
if exist('findPupilRayMex','file')==3
    findPupilHandle = @findPupilRayMex;
else
    findPupilHandle = @findPupilRay;
end

% Define a set of points on the perimeter of the stop
opts.Results.nStopPerimPoints = nStopPerimPoints;
opts.Results.stopPerimPhase = 0;
sg.eye = eye;
stopCoords = addStopPerimeter(sg,opts,[nan nan nan stopRadius]);

% Define the fieldRay for the passed parameters
fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);

% Detect if we are in a paraxial context. This occurs when the origin of
% the field ray is very close to the longitudinal axis. In this setting, we
% add a tiny bit angle to the field ray.
[SOR1,SOR2] = quadric.rayToAngles(quadric.coordsToRay([mean(stopCoords)',fieldRay(:,1)]));
if norm([SOR1,SOR2]) < paraxialThresh
    fieldRay = calcFieldRay([paraxialThresh/sqrt(2) paraxialThresh/sqrt(2)],rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);
end

% The location in object space is the first column of the field ray.
objectCoord = fieldRay(:,1);

% Create the optical systems
opticalSystemA = assembleOpticalSystem(eye,...
    'surfaceSetName','stopToMedium','cameraMedium',cameraMedium);
opticalSystemB = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToCamera','cameraMedium',cameraMedium);

% Assemble the static args for the findPupilRay
args = {convertEyeToWorldCoord(objectCoord'), ...
    eye.rotationCenters, ...
    opticalSystemA, ...
    opticalSystemB};

% Loop through the stop perimeter points. For each stop point, trace a ray
% from the perimeter to the object location.
p = nan(3,nStopPerimPoints);
u = nan(3,nStopPerimPoints);

for ii = 1:nStopPerimPoints
    
    % The (transposed) output ray for this stop perimeter location
    outputRay = ...
        findPupilHandle(stopCoords(ii,:), [0 0 0 nan], args{:})';
     
    % Store the ray components
    p(:,ii) = outputRay(:,1);
    u(:,ii) = outputRay(:,2);
    
end

% Some mechanics to shift the entrance window points back along their ray
% path to place the center of the entrance window on the longitudinal axis.
shifted = @(t) p - t.*u;
hvNorm = @(ep) norm(mean(ep(2:3,:),2));
myObj = @(t) hvNorm(shifted(t));

options = optimset('fminunc');
options.Display = 'off';
t = fminunc(myObj,0,options);

% These are our return variables
entranceWindowPerimeter = shifted(t);
entranceWindowCenter = mean(entranceWindowPerimeter,2);
radiusEntranceWindow = mean(vecnorm(entranceWindowPerimeter-entranceWindowCenter));

end