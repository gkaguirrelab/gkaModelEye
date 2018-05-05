function [outputRay, azimuths, elevations, intersectionCoords] = rayTraceCentered3DSurfaces(coordsInitial, angleInitial, opticalSystemIn, figureFlag)
% Returns the position and angle of a resultant ray w.r.t. the optical axis
%
% Syntax:
%  [outputRay, thetas, imageCoords, intersectionCoords] = rayTraceCenteredSurfaces(coordsInitial, thetaInitial, opticalSystemIn, figureFlag)
%
% Description:
%   This routine implements a 3D version of the generalized ray tracing
%   equations of:
%
%       Elagha, Hassan A. "Generalized formulas for ray-tracing and
%       longitudinal spherical aberration." JOSA A 34.3 (2017): 335-343.
%
%   The implementation assumes a set of elliptical surfaces, with each
%   surface having its center positioned on the optical axis. The initial
%   state of the ray is specified by its two-dimensional coordinates and by
%   the angle (theta) that it makes with the optical axis. By convention,
%   the optical axis is termed "z", and the orthogonal axis is termed
%   "height". Positive values of z are to the right. A theta of zero
%   indicates a ray that is parallel to the optical axis. Positive values
%   of theta correspond to the ray diverging to a position above the
%   optical axis. Each elliptical surface is specified by a center and a
%   radius in the z and h dimensions. The center must lie on the optical
%   axis; positive values place the center to the right of the origin of
%   the ray. A positive radius presents the ray with a convex surface; a
%   negative radius presents the ray with a concave surface. The output of
%   the routine is the position and angle at which the ray (or its reverse
%   projection) intersects the optical axis.
%
% Inputs:
%   coordsInitial         - A 2x1 orf 3x1 vector, with the values 
%                           corresponding to the z-position, height, and
%                           optionally depth of the initial position of the
%                           ray.
%   thetaInitial          - A scalar or 2x1 vector in radians that
%                           specifies the angle of the ray w.r.t. the
%                           optical axis in height plane and optionally in
%                           the depth plabne. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   opticalSystemIn       - An mx3, mx4, or mx5 matrix, where m is the 
%                           number of surfaces in the model, including the
%                           initial state of the ray. Each row contains the
%                           values:
%                               [center, radius, refractiveIndex]
%                           or
%                               [center, radiusZ, radiusH, refractiveIndex]
%                           or
%                               [center, radiusZ, radiusH, radiusK, refractiveIndex]
%                           that define an elliptical lens. The first row
%                           corresponds to the initial conditions of the
%                           ray. Thus, the refractive index value given in
%                           the first row specifies the index of the medium
%                           in which the ray arises. The center and radius
%                           values for the first row are ignored.
%   figureFlag            - Logical or structure. If logical, the true or
%                           false value controls whether the default style
%                           plot is created. If empty, it will be set to
%                           false. More control of the plot can be obtained
%                           by passing a structure with individual elements
%                           set to true or false.
%
% Outputs:
%   outputRay             - A 2x3 matrix that contains a unit vector that
%                           describes the location and vector direction of
%                           a ray that is the virtual image for this system
%                           for the input. The first columns contains the z
%                           positions, the second column the h(eight)
%                           position, and the third column the d(epth). The
%                           first row contains the coordinates of a point
%                           on the optic axis and the second row contains a
%                           point on a ray arising from the first point
%                           that has unit length. This vector originates
%                           from the point in space and has the same theta
%                           as the ray which emerges from the final surface
%                           for the input ray.
%   thetas                - A an mx2 ve in radians with the theta values at
%                           each surface.
%   imageCoords           - An mx3 matrix which provides at each surface
%                           the point at which the resultant ray (or its
%                           virtual extension) intersects the optical axis.
%                           The second column of this matrix will contain
%                           only zeros.
%   intersectionCoords    - An mx3 matrix that provides at each surface
%                           the point at which the ray intersects the
%                           surface.
%
% Examples:
%{
    %% Elagha 2017 numerical example
    % The paper provides a numerical example in section C which is
    % implemented here as an example. Compare the returned theta values
    % with those given on page 340, section C.
    clear figureFlag
    coords = [0 0];
    theta = deg2rad(17.309724);
    figureFlag=true;
    opticalSystem=[nan nan 1; 22 10 1.2; 9 -8 1; 34 12 1.5; 20 -10 1.0];
    [outputRay, thetas, phis, intersectionCoords] = rayTraceCentered3DSurfaces(coords, theta, opticalSystem, figureFlag);
    for ii=1:length(thetas)
        fprintf('theta%d: %f \n',ii-1,rad2deg(thetas(ii)));
    end
%}
%{
    %% Pupil through cornea
    % A model of the passage of a point on the pupil perimeter through
    % the axial cross-section of the cornea (units in mm)
    sceneGeometry = createSceneGeometry();
    [outputRay, thetas, phis, intersectionCoords] = rayTraceCentered3DSurfaces([sceneGeometry.eye.pupil.center(1) 2], [deg2rad(-15) 0], sceneGeometry.refraction.opticalSystem, true);
    % Compare the output to value calculated on April 26, 2018
    outputRayCached = [-0.052350606943506   0.882953893324114;   0.869847276998802   0.496235568373862];
    assert ( max(max(abs(outputRayCached - outputRay))) < 1e-6)
%}

%{
    %% Pupil through cornea, multiple points and rays
    sceneGeometry = createSceneGeometry();
    pupilRadius = 2;
    % Define FigureFlag as a structure, and set the new field to false so
    % that subsequent calls to the ray tracing routine will plot on the
    % same figure. Also, set the textLabels to false to reduce clutter
    figure
    clear figureFlag
    figureFlag.new = false;
    figureFlag.textLabels = false;
    for theta = -35:70:35
        for pupilRadius = -2:4:2
            rayTraceCenteredSurfaces([sceneGeometry.eye.pupil.center(1) pupilRadius], theta, sceneGeometry.refraction.opticalSystem.p1p2, figureFlag);
        end
    end
%}
%{
    %% Demo warnings
    coords = [0 0];
    opticalSystem=[nan nan 1.5; 20 10 1.0];
    % This ray will not intersect the surface. The function issues
    % warning and returns an empty outputRay
    theta = deg2rad(45);
    outputRay = rayTraceCenteredSurfaces(coords, theta, opticalSystem);
    % Make the index of refraction of the surface very high
    opticalSystem=[nan nan 5; 20 10 1.0];
    % This ray encounters total internal reflection. The function issues
    % warning and returns an empty outputRay
    theta = deg2rad(15);
    outputRay = rayTraceCenteredSurfaces(coords, theta, opticalSystem);
%}


%% Tell codegen to ignore these functions
% This is needed for compilation of the function as standalone mex file
coder.extrinsic('warning','legend','strseq')


%% Check the input
% If a value was not passed for figureFlag, set to false.
if nargin==3
    figureFlag.show = false;
    figureFlag.new = false;
    figureFlag.refLine = false;
    figureFlag.surfaces = false;
    figureFlag.imageLines = false;
    figureFlag.rayLines = false;
    figureFlag.finalUnitRay = true;
    figureFlag.textLabels = false;
    figureFlag.legend = false;
    figureFlag.zLim = [];
    figureFlag.hLim = [];
end

% A value was passed for figureFlag
if nargin==4
    % if the passed figureFlag is a structure, set all fields to true, and
    % then copy over the values of passed fields. This allows the use to
    % just specify values for some fields and causes the remainder to have
    % valid values.
    if isstruct(figureFlag)
        temp=figureFlag;
        clear figureFlag
        figureFlag.show = true;
        figureFlag.new = true;
        figureFlag.refLine = true;
        figureFlag.surfaces = true;
        figureFlag.imageLines = true;
        figureFlag.rayLines = true;
        figureFlag.finalUnitRay = true;
        figureFlag.textLabels = true;
        figureFlag.legend = true;
        figureFlag.zLim = [];
        figureFlag.hLim = [];
        names = fieldnames(temp);
        for nn = 1:length(names)
            figureFlag.(names{nn})=temp.(names{nn});
        end
    end
    % if figureFlag is logical, set all fields to true or false accordingly
    if islogical(figureFlag)
        if figureFlag
            clear figureFlag
            figureFlag.show = true;
            figureFlag.new = true;
            figureFlag.refLine = true;
            figureFlag.surfaces = true;
            figureFlag.imageLines = true;
            figureFlag.rayLines = true;
            figureFlag.finalUnitRay = true;
            figureFlag.textLabels = true;
            figureFlag.legend = true;
            figureFlag.zLim = [];
            figureFlag.hLim = [];
        else
            clear figureFlag
            figureFlag.show = false;
            figureFlag.new = false;
            figureFlag.refLine = false;
            figureFlag.surfaces = false;
            figureFlag.imageLines = false;
            figureFlag.rayLines = false;
            figureFlag.finalUnitRay = false;
            figureFlag.textLabels = false;
            figureFlag.legend = false;
            figureFlag.zLim = [];
            figureFlag.hLim = [];
        end
    end
end


%% Initialize variables and plotting
% outputRay set to empty
outputRay = [];
% strip the optical system of any rows which are all nans
opticalSystemIn=opticalSystemIn(sum(isnan(opticalSystemIn),2)~=size(opticalSystemIn,2),:);
% determine the number of surfaces and dimensionalty in which they are
% defined (i.e., spherical or elliptical)
nSurfaces = size(opticalSystemIn,1);
nDims = size(opticalSystemIn,2);

% Pre-allocate our loop variables; set the values for the first surface
% (initial position of ray)
aVals_p1p2 = zeros(nSurfaces,1);
aVals_p1p2(1,:) = 1;
aVals_p1p3 = zeros(nSurfaces,1);
aVals_p1p3(1,:) = 1;
curvature_p1p2 = zeros(nSurfaces,1);
curvature_p1p3 = zeros(nSurfaces,1);
curvatureCenters_p1p2 = zeros(nSurfaces,1);
curvatureCenters_p1p3 = zeros(nSurfaces,1);
intersectionCoords=zeros(nSurfaces,3);
if length(coordsInitial)==2
    intersectionCoords(1,:)=[coordsInitial 0];
else
    intersectionCoords(1,:)=coordsInitial;
end
azimuths = zeros(nSurfaces,1);
elevations = zeros(nSurfaces,1);
if length(angleInitial)==1
    azimuths(1)=angleInitial;
else
    azimuths(1)=angleInitial(1);
    elevations(1)=angleInitial(2);
end
relativeIndices = zeros(nSurfaces,1);
relativeIndices(1,:) = 1;

% Build the local optical system. Replace the center and radius of the
% first surface with the point of intersection of the initial ray with the
% optical axis, and set the radius to zero.
opticalSystem = zeros(nSurfaces,5);

% If the radii of each ellipsoidal surface is defined in one or two
% dimensions, copy the trailing value over to define the 3D surface. Thus,
% if a single radius is provided the system will model a sphere. If two
% radius values are provided the system will model an ellipsoid that is
% rotationally symmetric about the optical axis.
switch nDims
    case 3
        opticalSystem = [opticalSystemIn(:,1) opticalSystemIn(:,2) opticalSystemIn(:,2) opticalSystemIn(:,2) opticalSystemIn(:,3)];
    case 4
        opticalSystem = [opticalSystemIn(:,1) opticalSystemIn(:,2) opticalSystemIn(:,3) opticalSystemIn(:,3) opticalSystemIn(:,4)];
    case 5
        opticalSystem = opticalSystemIn;
end

% Initialize the figure
if figureFlag.show
    if figureFlag.new
        if figureFlag.legend
            figure
            subplot(3,3,1:6);
        else
            figure
        end
    else
        if figureFlag.legend
            subplot(3,3,1:6);
        else
        end
    end
    if nDims == 5
        subplot(3,3,4:6);
        hold on
        axis equal
        subplot(3,3,1:3);
    end
    hold on
    axis equal
    if ~isempty(figureFlag.zLim)
        xlim(figureFlag.zLim);
    end
    if ~isempty(figureFlag.hLim)
        ylim(figureFlag.hLim);
    end
end


%% Peform the ray trace
for ii = 2:nSurfaces

    % Obtain the coordinate at which the ray intersects the next surface,
    % and the radius of curvature of the surface at that point in the p1p2
    % and p1p3 planes
    [ intersectionCoords(ii,:), curvature_p1p2(ii), curvature_p1p3(ii), ellipseRadii_p1p2, ellipseRadii_p1p3 ] = rayIntersectEllipsoid( intersectionCoords(ii-1,:), azimuths(ii-1), elevations(ii-1), [opticalSystem(ii,2) opticalSystem(ii,3) opticalSystem(ii,4)], [opticalSystem(ii,1) 0 0] );
    % Check if the ray missed (or was tangenital to) the surface
    if isnan(curvature_p1p2(ii)) || isnan(curvature_p1p3(ii))
        warning('rayTraceCenteredSurfaces:nonIntersectingRay','Ray did not intersect surface %d. Returning.',ii);
        return
    end
    % Find the curvature center, which is the position along the optical
    % axis for this surface, given the curvature encountered by the ray.
    curvatureCenters_p1p2(ii) = opticalSystem(ii,1)-opticalSystem(ii,2)+curvature_p1p2(ii);
    curvatureCenters_p1p3(ii) = opticalSystem(ii,1)-opticalSystem(ii,3)+curvature_p1p3(ii);
    % The distance between the center of curvature of the current lens
    % surface and the center of curvature of the prior lens surface
    d_p1p2 = curvatureCenters_p1p2(ii)-curvatureCenters_p1p2(ii-1);
    d_p1p3 = curvatureCenters_p1p3(ii)-curvatureCenters_p1p3(ii-1);
    % The relative refractive index of the prior medium to the medium of
    % the surface that the ray is now impacting
    relativeIndices(ii)=opticalSystem(ii-1,end)/opticalSystem(ii,end);
    % Equation 54 of Elagha
    aVals_p1p2(ii) = ...
        (1/curvature_p1p2(ii))*(relativeIndices(ii-1).*aVals_p1p2(ii-1).*curvature_p1p2(ii-1)+d_p1p2.*sin(azimuths(ii-1)));
    aVals_p1p3(ii) = ...
        (1/curvature_p1p3(ii))*(relativeIndices(ii-1).*aVals_p1p3(ii-1).*curvature_p1p3(ii-1)+d_p1p3.*sin(elevations(ii-1)));
    % Check if the incidence angle is above the critical angle for the
    % relative refractive index at the surface interface.
    if abs((aVals_p1p2(ii)*relativeIndices(ii))) > 1 || abs((aVals_p1p3(ii)*relativeIndices(ii))) > 1
        warning('rayTraceCenteredSurfaces:criticalAngle','Angle of incidence for surface %d greater than critical angle. Returning.',ii);
        return
    end
    % Find the angle of the ray after it enters the current surface
    azimuths(ii) = azimuths(ii-1) - asin(aVals_p1p2(ii)) + asin(aVals_p1p2(ii).*relativeIndices(ii));
    elevations(ii) = elevations(ii-1) - asin(aVals_p1p3(ii)) + asin(aVals_p1p3(ii).*relativeIndices(ii));
    % Update the plot
    if figureFlag.show
        % add this lens surface
        if figureFlag.surfaces
            if nDims==5
                subplot(3,3,1:3);
                plotLensArc([opticalSystem(ii,1) ellipseRadii_p1p2])
                subplot(3,3,4:6);
                plotLensArc([opticalSystem(ii,1) ellipseRadii_p1p3])
            else
                plotLensArc(opticalSystem(ii,[1 2 3]))
            end
        end
        % plot the line for the path of the ray
        if figureFlag.rayLines
            if nDims==5
                subplot(3,3,1:3);
            plot([intersectionCoords(ii-1,1) intersectionCoords(ii,1)],[intersectionCoords(ii-1,2) intersectionCoords(ii,2)],'-r');
                subplot(3,3,4:6);
            plot([intersectionCoords(ii-1,1) intersectionCoords(ii,1)],[intersectionCoords(ii-1,3) intersectionCoords(ii,3)],'-r');
            else
            plot([intersectionCoords(ii-1,1) intersectionCoords(ii,1)],[intersectionCoords(ii-1,2) intersectionCoords(ii,2)],'-r');
            end
        end
    end
end


%% Finish and clean up
% Assemble an output which is the unit vector for the final ray
slope_theta = tan(azimuths(nSurfaces)+pi);
slope_phi = tan(elevations(nSurfaces)+pi);
norm_p1p2 = sqrt(slope_theta^2+1);
norm_p1p3 = sqrt(slope_theta^2+1);
outputRay = [intersectionCoords(nSurfaces,:); [intersectionCoords(nSurfaces,1)+(1/norm_p1p2) intersectionCoords(nSurfaces,2)+(slope_theta/norm_p1p2) intersectionCoords(nSurfaces,3)+(slope_phi/norm_p1p3)]];

% Complete the plot
if figureFlag.show
    % Plot the ray path, which is a triple-length output ray. Also, add a
    % mark to indicate the initial coordinates of the ray
    if figureFlag.rayLines
        if nDims==5
            subplot(3,3,1:3);
            finalRay = [intersectionCoords(nSurfaces,[1 2]); [intersectionCoords(nSurfaces,1)+(3/norm_p1p2) intersectionCoords(nSurfaces,2)+(3*slope_theta/norm_p1p2)]];
            plot([finalRay(1,1) finalRay(2,1)],[finalRay(1,2) finalRay(2,2)],'-r');
            plot(intersectionCoords(1,1),intersectionCoords(1,2),'xr');
            subplot(3,3,4:6);
            finalRay = [intersectionCoords(nSurfaces,[1 3]); [intersectionCoords(nSurfaces,1)+(3/norm_p1p3) intersectionCoords(nSurfaces,3)+(3*slope_phi/norm_p1p3)]];
            plot([finalRay(1,1) finalRay(2,1)],[finalRay(1,2) finalRay(2,2)],'-r');
            plot(intersectionCoords(1,1),intersectionCoords(1,3),'xr');
        else
            finalRay = [intersectionCoords(nSurfaces,[1 2]); [intersectionCoords(nSurfaces,1)+(3/norm_p1p2) intersectionCoords(nSurfaces,2)+(3*slope_theta/norm_p1p2)]];
            plot([finalRay(1,1) finalRay(2,1)],[finalRay(1,2) finalRay(2,2)],'-r');
            plot(intersectionCoords(1,1),intersectionCoords(1,2),'xr');
        end
    end
    % Plot the output unit ray vector
    if figureFlag.finalUnitRay
        if nDims==5
            subplot(3,3,1:3);
            plot([outputRay(1,1) outputRay(2,1)],[outputRay(1,2) outputRay(2,2)],'-g');
            subplot(3,3,4:6);
            plot([outputRay(1,1) outputRay(2,1)],[outputRay(1,3) outputRay(2,3)],'-g');
        else
            plot([outputRay(1,1) outputRay(2,1)],[outputRay(1,2) outputRay(2,2)],'-g');
        end
    end
    % Replot the refline
    if figureFlag.refLine
        yl = [min([intersectionCoords(:,1); opticalSystem(:,1)]),max([intersectionCoords(:,1); opticalSystem(:,1)])];
        yl(1) = yl(1)-(yl(1)/3)*sign(yl(1));
        yl(2) = yl(2)+yl(1)/3*sign(yl(1));
        if nDims==5
            subplot(3,3,1:3);
            plot(yl,[0,0],'-k');
            subplot(3,3,4:6);
            plot(yl,[0,0],'-k');
        else
            plot(yl,[0,0],'-k');
        end
    end
    % Add some labels
    if figureFlag.textLabels
        if nDims==5
            subplot(3,3,1:3);
            yl = max(abs(intersectionCoords(:,2)));
            yPos = zeros(nSurfaces,1)-(yl/2);
            plot(opticalSystem(:,1),zeros(nSurfaces,1),'+k');
            labels = strseq('c',1:1:nSurfaces);
            text(opticalSystem(:,1),yPos,labels,'HorizontalAlignment','center');
            subplot(3,3,4:6);
        end
        yl = max(abs(intersectionCoords(:,2)));
        yPos = zeros(nSurfaces,1)-(yl/2);
        plot(opticalSystem(:,1),zeros(nSurfaces,1),'+k');
        labels = strseq('c',1:1:nSurfaces);
        text(opticalSystem(:,1),yPos,labels,'HorizontalAlignment','center');
    end
    % Add a legend
    if figureFlag.legend
        hSub = subplot(3,3,8);
        plot(nan, nan,'-r');
        hold on
        plot(nan, nan,'-b');
        plot(nan, nan,'-g');
        set(hSub, 'Visible', 'off');
        legend({'ray path','final virtual ray','output unit ray vector'},'Location','north', 'Orientation','vertical');
    end
    hold off
end

end % function


%% LOCAL FUNCTIONS

function plotLensArc(opticalSystem)
% Local function to handle plotting lens surfaces
ang=pi/2:0.01:3*pi/2;
xp=opticalSystem(2)*cos(ang);
yp=opticalSystem(3)*sin(ang);
plot(opticalSystem(1)+xp,yp,'-k');
end


function [ coordsOut, curvature_p1p2, curvature_p1p3, ellipseRadii_p1p2, ellipseRadii_p1p3 ] = rayIntersectEllipsoid( coordsIn, azimuth, elevation, ellipsoidRadii, ellipsoidCenter )
% Returns coords and curvature at intersection of a ray and an ellipsoid
%
% Syntax:
%  [intersectionCoords, curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
%
% Description:
%   Implements trigonometric operations to identify the point at which a
%   line intersects an ellipse, and returns the curvature of the ellipse at
%   the point of contact.
%
% Inputs:
%   coordsInitial         - 2x1 vector, with the values corresponding to
%                           the z-position and height of the initial
%                           position of the ray.
%   theta                 - Scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   ellipseCenterZ        - Scalar in mm. The position of the center of the
%                           elliptical surface on the z-axis
%   ellipseRadii          - A 2x1 vector in mm, with the values
%                           corresponding to the radii of the ellipitical
%                           surface along the z and height axes,
%                           respectively.
%
% Outputs:
%   intersectionCoords    - 2x1 vector with the values corresponding to the
%                           z-position and height at which the ray strikes
%                           the elliptical surface
%   curvature             - Scalar. The radius of curvature of the
%                           elliptical surface at the point of intersection
%                           with the ray.
%
% Examples:
%{
    [intersectCoords,curvature] = calcEllipseIntersect([-3.7,0], 0, -9.1412, [9.1412, 8.4277] )
%}


% Initialize return variables
coordsOut = [nan, nan, nan];
curvature_p1p2 = nan;
curvature_p1p3 = nan;
ellipseRadii_p1p2 = [nan nan];
ellipseRadii_p1p3 = [nan nan];

% Store the sign of the radius values. They radii must have the same sign.
radiiSign = sign(ellipsoidRadii(1));
if ~all(radiiSign == radiiSign(1))
    error('rayTraceCenteredSurfaces:incompatibleConvexity','The radii of the elliptical lens surface must have the same sign.');
end
% Convert the radii to their absolute values
ellipsoidRadii = abs(ellipsoidRadii);

% Define a ray as P + tu
% Where P is the point of origin of the ray, u is the unit vector direction
% of the ray, and t is the weight on that unit vector direction

% convert aziuth and elevation to theta and phi
theta=acos(cos(elevation)*cos(azimuth))*sign(azimuth);
phi=atan(tan(elevation)/sin(azimuth));

ray = createLine3d(coordsIn([1 3 2]),theta,phi);
ray = ray([1 3 2 6 4 5]);
unitSphere = [0 0 0 1];

% Consider the transformation of the unit sphere to an ellipsoid by:
%  - Scaling in the p1, p2, and p3 directions
%  - Rotation about the p1, p2, and p3 directions (not used here)
%  - Translation (we are only using shifts along the optical axis)
% We combine these in the 4x4 transformation matrix K.
% We then find the intersection of the ray inv(K)P + t*inv(K)u with the unit sphere.

% Construct K = translate * rotate * scale
scale = [ellipsoidRadii(1) 0 0 0; 0 ellipsoidRadii(2) 0 0; 0 0 ellipsoidRadii(3) 0; 0 0 0 1];
rotate = eye(4,4);
translate = eye(4,4);
translate(1:3,4) = ellipsoidCenter;
K = translate * rotate * scale;
invK = K\eye(4,4);

% transform the ray
kP = invK*[ray(1:3) 1]';
ku = invK*[ray(4:6) 0]';
kRay = [kP(1:3)' ku(1:3)'];

% Find the transformed intersection points on the unit sphere
intersectionPoints = intersectLineSphere(kRay, unitSphere);

% Detect if we have a tangential ray
if sum(abs(intersectionPoints(1,:)-intersectionPoints(2,:)))==0
    return
end

% Detect if we have a non-intersecting ray
if isnan(intersectionPoints(1,1))
    return
end

% transform the intersection points back to the original coordinate frame
coordsOutA = K*[intersectionPoints(1,:) 1]';
coordsOutB = K*[intersectionPoints(2,:) 1]';
coordsOutAB = [coordsOutA(1:3)'; coordsOutB(1:3)'];

% If the radiiSign is positive, report the coordinates on the left-hand
% side of the ellipse, otherwise report the coordinates on the right
if radiiSign<0
    [~,rightIdx] = max(coordsOutAB(:,1));
    coordsOut = coordsOutAB(rightIdx,:);
else
    [~,leftIdx] = min(coordsOutAB(:,1));
    coordsOut = coordsOutAB(leftIdx,:);
end

% Now find the curvature of the ellipsoid at the intersection. To do so, we
% obtain the radii of the ellipse in the orthogonal planes of the rotated
% ray for the theta (p1p2) and phi (p1p3) angles.


% The radii of the ellipse that lies within the p1p2 plane when the ray is
% diverging from the optical axis into the p3 dimension by angle elevation
A = tan(elevation); B = 0; C = 1;
D = coordsOut(2)/((coordsOut(1)-ellipsoidCenter(1))*A);
if isinf(D) || isnan(D)
    D = 0;
end
ellipseRadii_p1p2 = [complex(0) complex(0)];
[ellipseRadii_p1p2(1),ellipseRadii_p1p2(2)]=EllipsoidPlaneIntersection(A,B,C,0,ellipsoidRadii(1),ellipsoidRadii(2),ellipsoidRadii(3));
ellipseRadii_p1p2 = real(ellipseRadii_p1p2);

% The radii of the ellipse that lies within the p1p3 plane when the ray is
% diverging from the optical axis into the p2 dimension by angle azimuth
A = tan(azimuth); B = 1; C = 0;
D = coordsOut(3)/((coordsOut(1)-ellipsoidCenter(1))*A);
if isinf(D) || isnan(D)
    D = 0;
end
ellipseRadii_p1p3 = [complex(0) complex(0)];
[ellipseRadii_p1p3(1),ellipseRadii_p1p3(2)]=EllipsoidPlaneIntersection(A,B,C,0,ellipsoidRadii(1),ellipsoidRadii(2),ellipsoidRadii(3));
ellipseRadii_p1p3 = real(ellipseRadii_p1p3);

% Calculate the radius of curvature encountered by the ray.
% p1p2 plane
t = real(acos(complex((coordsOut(1)-ellipsoidCenter(1))/ellipseRadii_p1p2(1))));
curvature_p1p2 = radiiSign*((ellipseRadii_p1p2(1)^2*sin(t)^2 + ellipseRadii_p1p2(2)^2*cos(t)^2)^(3/2))/(ellipseRadii_p1p2(1)*ellipseRadii_p1p2(2));
% p1p3 plane
t = real(acos(complex((coordsOut(1)-ellipsoidCenter(1))/ellipseRadii_p1p3(1))));
curvature_p1p3 = radiiSign*((ellipseRadii_p1p3(1)^2*sin(t)^2 + ellipseRadii_p1p3(2)^2*cos(t)^2)^(3/2))/(ellipseRadii_p1p3(1)*ellipseRadii_p1p3(2));

% Adjust the ellipse radii for the sign of the input radii and return these
ellipseRadii_p1p2 = ellipseRadii_p1p2 * radiiSign;
ellipseRadii_p1p3 = ellipseRadii_p1p3 * radiiSign;


end

