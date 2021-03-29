function [incidentNode,emergentNode,incidentRays,emergentRays] = calcNodes(opticalSystem,longitudinalFieldAngle,rayOriginDistance,bundleCount,cameraMedium)
% The approximate incident and emergent nodal points of an eye
%
% Syntax:
%  [incidentNode,emergentNode,incidentRays,emergentRays] = calcNodes(opticalSystem,longitudinalFieldAngle,rayOriginDistance,bundleCount,cameraMedium)
%
% Description
%   For the aspheric, astigmatic optical system of the eye, single nodal
%   points do not exist. Nonetheless, an approximation of the nodal points
%   is useful. Given a model eye, this routine finds nodal rays from
%   different points in the optical field. The bundle is then examined to
%   identify the axial position at which the cross-sectional area of the
%   incident and emergent ray bundles is smallest (i.e., the "waist" of the
%   bundle). These are the incident and emergent nodal ellipses. The
%   intersection of the plane of these ellipses with the optical axis of
%   the eye is taken as the location of the incident and emergent nodal
%   points.
%
%   These ideas are discussed in:
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
% Inputs:
%   opticalSystem         - Either an eye structure (from which a
%                           "mediumToRetina" optical system in air will be
%                           derived), or an mx19 matrix, where m is set by
%                           the key value opticalSystemNumRows. Each row
%                           contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%   longitudinalFieldAngle - Scalar. The optical field angle that will be
%                           used to define the bundle of rays.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the rays from the longitudinal axis origin.
%   bundleCount           - Scalar. The number of rays to examine.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air' if not provided.
%
% Outputs:
%   incidentNode          - 3x1 vector which provides the coordinates in
%                           the [p1 p2 p3] dimensions of the approximate
%                           incident nodal point of the eye.
%   emergentNode          - 3x1 vector which provides the coordinates in
%                           the [p1 p2 p3] dimensions of the approximate
%                           emergent nodal point of the eye.
%   incidentRays          - A cell array, with each cell containing a 3x2
%                           matrix that provides the incident nodal rays.
%   emergentRays          - A cell array, with each cell containing a 3x2
%                           matrix that provides the emergent nodal rays.
%
% Examples:
%{
    % Find the approximate nodes of a model eye
    eye = modelEyeParameters();
    [incidentNode,emergentNode] = calcNodes(eye);
    fprintf('The longitudinal positions of the approximate incident and emergent nodes are %0.2f and %0.2f mm \n',incidentNode(1),emergentNode(1));
%}
%{
    % Find and display the nodes on the model eye
    eye = modelEyeParameters();
    % Obtain the optical center and ray bundle
    [incidentNode,~,incidentRays,~] = calcNodes(eye);
    % Scale up the incidentRays
    extendedRays = cellfun(@(x) [x(:,1),x(:,1)+x(:,2)*1520],incidentRays,'UniformOutput',false);
    % Plot the optical system
    plotOpticalSystem(eye,'surfaceAlpha', 0.05);
    % Add the rays to the plot
    for ii=1:length(incidentRays)
        plotOpticalSystem('newFigure',false,'rayPath',extendedRays{ii});
    end
    xlim([-25 10]);
    plot3(incidentNode(1),incidentNode(2),incidentNode(3),'*k');
%}


arguments
    opticalSystem {mustBeOpticalSystemCapable}
    longitudinalFieldAngle (1,1) {mustBeNumeric} = 14.3; % Following Harris Example 3
    rayOriginDistance (1,1)  {mustBeNumeric} = 1500
    bundleCount (1,1)  {mustBeNumeric} = 16
    cameraMedium = 'air'
end

% Create the optical system
opticalSystem = parseOpticalSystemArgument(opticalSystem,'mediumToRetina',cameraMedium);

% Find the optical axis for this system
opticalAxisRay = calcOpticalAxis(opticalSystem, rayOriginDistance);

% Prepare variables for the loop
rayPaths = cell(1,bundleCount);

% Loop over longitudinalFieldAngles and send off the ray bundle
for bb = 1:bundleCount
    opticalFieldOrigin = [ ...
        sind(bb*360/bundleCount)*longitudinalFieldAngle, ...
        cosd(bb*360/bundleCount)*longitudinalFieldAngle, ...
        ];
    rayPaths{bb} = calcNodalRayFromField(opticalSystem,opticalFieldOrigin,rayOriginDistance);
end

% Obtain the sets of incident and emergent rays
incidentRays = cellfun(@(x) quadric.coordsToRay(x(:,1:2)),rayPaths,'UniformOutput',false);
emergentRays = cellfun(@(x) quadric.coordsToRay(x(:,end-1:end)),rayPaths,'UniformOutput',false);

% If we encountered bad ray traces, there will be nans in the emergent
% rays. If so, exit with a warning
if any(cellfun(@(x) any(isnan(x(:))),emergentRays))
    warning('calcNodes:badRayTrace','Unable to conduct ray trace. Reduce longitudinalFieldAngle, or check for invalid opticalSystem');
    return
end

% Obtain the set of unique pairs of the rays
[xSub,ySub]=ind2sub([bundleCount bundleCount],find(triu(ones(bundleCount,bundleCount),1)));

% We find the closest point between each unique pair of incident and
% emergent rays
incidentCenters = [];
emergentCenters = [];
for ss=1:length(xSub)
    incidentCenters(end+1,:) = quadric.distanceRays(incidentRays{xSub(ss)},incidentRays{ySub(ss)});
    emergentCenters(end+1,:) = quadric.distanceRays(emergentRays{xSub(ss)},emergentRays{ySub(ss)});
end

% Find the linear fit of a plane to the set of incident centers
incidentMeanCenter = mean(incidentCenters,1);
R = incidentCenters-incidentMeanCenter;
V = eig(R'*R);
n = V(:,1);

% Find the intersection of this incident nodal plane with the optical axis
incidentNode = quadric.intersectRayPlane(n,incidentMeanCenter',opticalAxisRay);

% Find the linear fit of a plane to the set of emergent centers
emergentMeanCenter = mean(emergentCenters,1);
R = emergentCenters-emergentMeanCenter;
V = eig(R'*R);
n = V(:,1);

% Find the intersection of this incident nodal plane with the optical axis
emergentNode = quadric.intersectRayPlane(n,emergentMeanCenter',opticalAxisRay);


end
