function [nodalPointCoord, outputRays, rayPaths] = calcEffectiveNodalPoint(eye,cameraMedium)
% Returns the effective nodal point for an eye in a given medium
%
% Syntax:
%  [nodalPointCoord, outputRays, rayPaths] = calcEffectiveNodalPoint(eye,cameraMedium)
%
% Description
%   The nodal points of a lens have the property that a ray aimed at the
%   first point will then emerge from second point with the same angle
%   relative to the optical axis. For the aspheric, astigmatic optical
%   system of the eye, single nodal points do not exist. Nonetheless, an
%   approximation of the nodal point is useful. Given a model eye, this
%   routine examines a bundle of rays arising from different points in the
%   retina. For each point, the ray is found that exits the corneal surface
%   at the same angles (relative to the optical axis) with which it arose
%   from the retina. The ray bundle is then examined to identify the axial
%   position at which the cross-sectional area of the ray bundle is
%   smallest (i.e., the "waist" fo the bundle). The center of the waist at
%   this location is returned as the effective nodal point coordinate.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air' if not provided.
%
% Outputs:
%   nodalPointCoord       - 3x1 vector which provides the coordinates in
%                           the [p1 p2 p3] dimensions of the effective
%                           nodal point center.
%   outputRays            - A cell array, with each cell containing a
%                           3x2 matrix that specifies a ray as a unit 
%                           vector.
%   rayPaths              - A cell array, with each cell containing a 3xm
%                           matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.                     
%
% Examples:
%{
    % Find the nodal point of a model eye
    eye = modelEyeParameters();
    nodalPointCoord = calcEffectiveNodalPoint(eye,'air')
%}
%{
    % Find and display the effective nodal point on the model eye
    sceneGeometry = createSceneGeometry();
    % Plot the optical system
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'surfaceAlpha', 0.05);    
    % Obtain the effective nodal point and ray bundle
    [nodalPointCoord, outputRays, rayPaths] = calcEffectiveNodalPoint(sceneGeometry.eye);
    % Add the rays to the plot
    for ii=1:length(outputRays)
        plotOpticalSystem('newFigure',false,'outputRay',outputRays{ii},'rayPath',rayPaths{ii});
    end
    plot3(nodalPointCoord(1),nodalPointCoord(2),nodalPointCoord(3),'*k');
%}

% Parse the inputs
if nargin<1
    error('Invalid number of input arguments');
end

if nargin==1
    cameraMedium = 'air';
end

% Define the output variable
nodalPointCoord = zeros(3,1);

% Assemble the optical system
opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName','retinaToCamera', 'cameraMedium', cameraMedium );

% Predefine some variables for use in the upcoming loop
outputRays={};
rayPaths={};

% Define some options for the fmincon call in the loop
options = optimoptions(@fmincon,...
    'Display','off');

% Loop over a few locations in the retinal surface, specified in
% ellipsoidal coordinates
for beta = [-85,-75]
    for omega = -180:40:180
        % Get this retinal coordinate
        coord = quadric.ellipsoidalGeoToCart( [beta, omega, 0], eye.retina.S );
        % Define an error function that reflects the difference in angles
        % from the initial ray and the output ray from the optical system
        myError = @(x) calcOffsetFromParallel(opticalSystem,assembleInputRay(coord,x(1),x(2)));
        % Supply an x0 guess as the ray that connects the retinal point
        % with the center of the pupil
        [~, angle_p1p2, angle_p1p3] = quadric.angleRays( [0 0 0; 1 0 0]', quadric.normalizeRay([coord'; eye.pupil.center-coord']') );
        angle_p1p2 = deg2rad(angle_p1p2);
        angle_p1p3 = -deg2rad(angle_p1p3);
        % Perform the search
        inputRayAngles = fmincon(myError,[angle_p1p2 angle_p1p3],[],[],[],[],[-pi/2,-pi/2],[pi/2,pi/2],[],options);
        % Calculate and save the outputRay and the raypath
        [outputRays{end+1},rayPaths{end+1}] = rayTraceQuadrics(assembleInputRay(coord,inputRayAngles(1),inputRayAngles(2)), opticalSystem);
    end
end

% Remove any NaNs from the rayPaths. This happens when a given ray misses
% one of the shells of the gradient lens model. The presence of a nan in
% the rayPath vector is expected and proper, but the nan entry must be
% removed to allow the fminsearch operation below to proceed.
nanFreeRayPaths = cellfun(@(x) x(:,~isnan(x(1,:))),rayPaths,'UniformOutput',false);

% Find the waist of the ray bundle along the axial (p1) dimension
bundleArea = @(p1) range(cellfun(@(x) interp1(x(1,:),x(3,:),p1,'linear'),nanFreeRayPaths)) * ...
    range(cellfun(@(x) interp1(x(1,:),x(2,:),p1,'linear'),nanFreeRayPaths));
nodalPointCoord(1) = fminsearch(bundleArea,-6);

% Find the geometric center of the waist in the p2 and p3 dimensions
nodalPointCoord(2) = mean(cellfun(@(x) interp1(x(1,:),x(2,:),nodalPointCoord(1),'linear'),nanFreeRayPaths));
nodalPointCoord(3) = mean(cellfun(@(x) interp1(x(1,:),x(3,:),nodalPointCoord(1),'linear'),nanFreeRayPaths));


end

% Local function. Converts angles relative to the optical axis to a unit
% vector ray.
function inputRay = assembleInputRay(p,angle_p1p2,angle_p1p3)
u = [1; tan(angle_p1p2); tan(angle_p1p3)];
u = u./sqrt(sum(u.^2));
inputRay = [p, u];
end

% Local function. Performs the ray trace through the optical system of the
% eye and then calculates the angle between the initial ray and the output
% ray.
function angleError = calcOffsetFromParallel(opticalSystem,inputRay)
    exitRay = rayTraceQuadrics(inputRay, opticalSystem);
    angleError = quadric.angleRays( inputRay, exitRay );
end