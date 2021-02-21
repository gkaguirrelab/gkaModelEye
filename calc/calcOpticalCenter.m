function [opticalCenter,nodes,rayPaths,incidentNodes,emergentNodes] = calcOpticalCenter(eye,cameraMedium)
% Returns the optical center and nodes for an eye in a given medium
%
% Syntax:
%  [opticalCenterCoord, outputRays, rayPaths] = calcOpticalCenter(eye,cameraMedium)
%
% Description
%   The nodal points of a lens have the property that a ray aimed at the
%   first point will emerge from second point with the same angle relative
%   to the optical axis. For the aspheric, astigmatic optical system of the
%   eye, single nodal points do not exist. Nonetheless, an approximation of
%   the nodal points is useful. Given a model eye, this routine finds the
%   nodal ray from different points in the visual field. The bundle is then
%   examined to identify the axial position at which the cross-sectional
%   area of the ray bundle is smallest (i.e., the "waist" of the bundle).
%   The center of the waist at this location is returned as the effective
%   nodal point coordinate, or more properly termed, the optical center.
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
    opticalCenter = calcOpticalCenter(eye,'air');
    fprintf('The optical center is %0.2f mm axial distance from the cornea, and %0.2f mm from the retinal vertex\n',opticalCenter(1),opticalCenter(1)-eye.landmarks.vertex.coords(1));
%}
%{
    % Find and display the optical center on the model eye
    sceneGeometry = createSceneGeometry();
    % Obtain the optical center and ray bundle
    [opticalCenter,nodes,rayPaths] = calcOpticalCenter(sceneGeometry.eye);
    % Plot the optical system
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.mediumToRetina,'addLighting',true,'surfaceAlpha', 0.05);    
    % Add the rays to the plot
    for ii=1:length(rayPaths)
        plotOpticalSystem('newFigure',false,'rayPath',rayPaths{ii});
    end
    xlim([-25 10]);
    plot3(opticalCenter(1),opticalCenter(2),opticalCenter(3),'*k');
%}


% Parse the inputs
if nargin<1
    error('Invalid number of input arguments');
end

if nargin==1
    cameraMedium = 'air';
end

% Check if we have a compiled version of findNodalRay
if exist('findNodalRayMex','file')==3
    findNodeHandle = @findNodalRayMex;
else
    findNodeHandle = @findNodalRay;
end

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Prepare variables for the loop
rayPaths = {};
incidentNodes = [];
emergentNodes = [];
% Loop over locations in the visual field, defined by angle (w.r.t. the
% un-rotated corneal apex for convenience) and distance
for rayOriginDistance = 100:250:1500
    for horiz = -60:15:50
        for vert = -60:15:50
            % The calculation doesn't work for a ray along the optical axis
            if horiz==0 && vert == 0
                continue
            end
            rayOrigin = quadric.anglesToRay([0;0;0],horiz,vert).*rayOriginDistance;
            rayOrigin = rayOrigin(:,2);
            [rayPath,nodalPoints,errors] = findNodeHandle(rayOrigin',opticalSystem);
            if errors(1)<1e-3
                rayPaths{end+1}=rayPath;
                incidentNodes(:,end+1)=nodalPoints(:,1);
                emergentNodes(:,end+1)=nodalPoints(:,2);
            end
        end
    end
end


% Remove any NaNs from the rayPaths. This happens when a given ray misses
% one of the shells of the gradient lens model. The presence of a nan in
% the rayPath vector is expected and proper, but the nan entry must be
% removed to allow the fminsearch operation below to proceed.
nanFreeRayPaths = cellfun(@(x) x(:,~isnan(x(1,:))),rayPaths,'UniformOutput',false);

% In some cases the optical system will have a ray that begins at the first
% surface in the system, in which case the first two positions in the ray
% path will be identical. Remove these to avoid upsetting the interp
% operation that follows.
samePosTol = 1e-4;
uniquePositionRayPaths = cellfun(@(x) x(:,logical([sum((diff(x')').^2)>samePosTol,1])),nanFreeRayPaths,'UniformOutput',false);

% Find the waist of the ray bundle along the axial (p1) dimension
bundleArea = @(p1) range(cellfun(@(x) interp1(x(1,:),x(3,:),p1,'linear'),uniquePositionRayPaths)) * ...
    range(cellfun(@(x) interp1(x(1,:),x(2,:),p1,'linear'),uniquePositionRayPaths));
opticalCenter(1) = fminsearch(bundleArea,-6);

% Find the geometric center of the waist in the p2 and p3 dimensions
opticalCenter(2) = mean(cellfun(@(x) interp1(x(1,:),x(2,:),opticalCenter(1),'linear'),uniquePositionRayPaths));
opticalCenter(3) = mean(cellfun(@(x) interp1(x(1,:),x(3,:),opticalCenter(1),'linear'),uniquePositionRayPaths));

% Find the median position of the nodes
nodes = [median(incidentNodes,2) median(emergentNodes,2)];

end
