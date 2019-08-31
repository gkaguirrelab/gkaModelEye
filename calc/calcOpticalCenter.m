function [opticalCenterCoord, outputRays, rayPaths] = calcOpticalCenter(eye,cameraMedium)
% Returns the optical center for an eye in a given medium
%
% Syntax:
%  [opticalCenterCoord, outputRays, rayPaths] = calcOpticalCenter(eye,cameraMedium)
%
% Description
%   The nodal points of a lens have the property that a ray aimed at the
%   first point will then emerge from second point with the same angle
%   relative to the optical axis. For the aspheric, astigmatic optical
%   system of the eye, single nodal points do not exist. Nonetheless, an
%   approximation of the nodal points is useful. Given a model eye, this
%   routine examines a bundle of rays arising from different points in the
%   retina. For each point, the ray is found that exits the corneal surface
%   at the same angles (relative to the optical axis) with which it arose
%   from the retina (a "nodal ray"). The ray bundle is then examined to
%   identify the axial position at which the cross-sectional area of the
%   ray bundle is smallest (i.e., the "waist" fo the bundle). The center of
%   the waist at this location is returned as the effective nodal point
%   coordinate, or more properly termed, the optical center.
%
%   These ideas are discussed in:
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air' if not provided.
%
% Outputs:
%   opticalCenterCoord    - 3x1 vector which provides the coordinates in
%                           the [p1 p2 p3] dimensions of the approximate
%                           optical center of the eye.
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
    % Find the optical center of a model eye
    eye = modelEyeParameters();
    opticalCenterCoord = calcOpticalCenter(eye,'air');
    fprintf('The optical center is %0.2f mm axial distance from the cornea, and %0.2f mm from the retinal vertex\n',opticalCenterCoord(1),opticalCenterCoord(1)-eye.landmarks.vertex.coords(1));
%}
%{
    % Find and display the optical center on the model eye
    sceneGeometry = createSceneGeometry();
    % Plot the optical system
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'surfaceAlpha', 0.05);    
    % Obtain the optical center and ray bundle
    [opticalCenterCoord, outputRays, rayPaths] = calcOpticalCenter(sceneGeometry.eye);
    % Add the rays to the plot
    for ii=1:length(outputRays)
        plotOpticalSystem('newFigure',false,'outputRay',outputRays{ii},'rayPath',rayPaths{ii});
    end
    plot3(opticalCenterCoord(1),opticalCenterCoord(2),opticalCenterCoord(3),'*k');
%}

% Parse the inputs
if nargin<1
    error('Invalid number of input arguments');
end

if nargin==1
    cameraMedium = 'air';
end

% Define the output variable
opticalCenterCoord = zeros(3,1);

% Predefine some variables for use in the upcoming loop
outputRays={};
rayPaths={};

% Loop over a few locations on the retinal surface, specified in
% ellipsoidal coordinates
for beta = [-85,-75]
    for omega = -180:40:180
        % Get this retinal coordinate
        X = quadric.ellipsoidalGeoToCart( [beta, omega, 0], eye.retina.S );
        % Get the nodal ray
        [outputRays{end+1},rayPaths{end+1}] = calcNodalRay(eye,[],X,cameraMedium);
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
opticalCenterCoord(1) = fminsearch(bundleArea,-6);

% Find the geometric center of the waist in the p2 and p3 dimensions
opticalCenterCoord(2) = mean(cellfun(@(x) interp1(x(1,:),x(2,:),opticalCenterCoord(1),'linear'),nanFreeRayPaths));
opticalCenterCoord(3) = mean(cellfun(@(x) interp1(x(1,:),x(3,:),opticalCenterCoord(1),'linear'),nanFreeRayPaths));


end
