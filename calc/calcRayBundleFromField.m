function [retinaPoints,rayPaths] = calcRayBundleFromField(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,nStopPerimPoints,paraxialThresh,cameraMedium)
% Directs a bundle of rays at the entrance window of the eye
%
% Syntax:
%  [retinaPoints,rayPaths] = calcRayBundleFromField(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,nStopPerimPoints,paraxialThresh,cameraMedium)
%
% Description
%
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
%   retinaPoints          - 3 x nStopPerimPoints. The set of points on the
%                           retinal surface
%   rayPaths              - Cell array of length nStopPerimPoints. Each
%                           cell holds a ray path.
%
% Examples:
%{
    % Define a default model eye
    accommodation = 1;
    stopRadius = 2;
    eye = modelEyeParameters('spectralDomain','vis','accommodation',accommodation);
    % Pick a field location
    rayOriginDistance = 1000/accommodation;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    distanceReferenceCoord = calcPrincipalPoint(eye);
    fieldAngularPosition = [0 0];
    % Trace the bundle
    [retinaPoints,rayPaths] = calcRayBundleFromField(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);
    % Plot the optical system
    opticalSystem = assembleOpticalSystem(eye,'surfaceSetName','mediumToRetina');
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'surfaceAlpha', 0.05);
    % Add the rays to the plot
    for ii=1:length(rayPaths)
        plotOpticalSystem('newFigure',false,'rayPath',rayPaths{ii});
    end
    xlim([-25 10]);
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0,0]
    stopRadius (1,1) {mustBeNumeric} = 2
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    nStopPerimPoints (1,1) {mustBeNumeric} = 16
    paraxialThresh (1,1) {mustBeNumeric} = 1e-3
    cameraMedium = 'air'
end


% Create the optical system
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Find the entrance window
[entranceWindowCenter,objectCoord,entranceWindowPerimeter] = ...
    calcEntranceWindow(eye,fieldAngularPosition,stopRadius,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,nStopPerimPoints,paraxialThresh,cameraMedium);

% Direct rays from the field point to the borders of the entrance window
retinaPoints = [];
rayPaths = {};

% Fill the window by scaling the set of rays towards the window center
for ss = [0,0.25,0.5,1.0]
    for ii = 1:size(entranceWindowPerimeter,2)
        windowTarget = entranceWindowPerimeter(:,ii);
        windowTarget = entranceWindowCenter + ss.*(windowTarget - entranceWindowCenter);
        if ss==0 && ii~=1
            continue
        end
        bundleRay = quadric.coordsToRay([objectCoord,windowTarget]);
        [outputRay, rayPaths{end+1}] = rayTraceQuadrics(bundleRay, opticalSystem);
        retinaPoints(:,end+1) = outputRay(:,1);
    end
end

end