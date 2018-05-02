function [outputRay, thetas, imageCoords, intersectionCoords] = rayTraceCenteredSurfaces(coordsInitial, thetaInitial, opticalSystemIn, figureFlag)
% Returns the position and angle of a resultant ray w.r.t. the optical axis
%
% Syntax:
%  [outputRay, thetas, imageCoords, intersectionCoords] = rayTraceCenteredSurfaces(coordsInitial, thetaInitial, opticalSystemIn, figureFlag)
%
% Description:
%   This routine implements the 2D generalized ray tracing equations of:
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
%   coordsInitial         - A 2x1 vector, with the values corresponding to
%                           the z-position and height of the initial
%                           position of the ray.
%   thetaInitial          - A scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   opticalSystemIn       - An mx3 or mx4 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           state of the ray. Each row contains the values:
%                               [center, radius, refractiveIndex]
%                           or
%                               [center, radiusZ, radiusH, refractiveIndex]
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
%   outputRay             - A 2x2 matrix that contains a unit vector that
%                           describes the location and vector direction of
%                           a ray that is the virtual image for this system
%                           for the input. The first columns contains the z
%                           positions and the second column the h(eight)
%                           position. The first row contains the
%                           coordinates of a point on the optic axis and
%                           the second row contains a point on a ray
%                           arising from the first point that has unit
%                           length. This vector originates from the point
%                           in space and has the same theta as the ray
%                           which emerges from the final surface for the
%                           input ray.
%   thetas                - A scalar in radians
%   imageCoords           - An mx2 matrix which provides at each surface
%                           the point at which the resultant ray (or its
%                           virtual extension) intersects the optical axis.
%                           The second column of this matrix will contain
%                           only zeros.
%   intersectionCoords    - An mx2 matrix which provides at each surface
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
    [outputRay, thetas, imageCoords] = rayTraceCenteredSurfaces(coords, theta, opticalSystem, figureFlag);
    for ii=1:length(thetas)
        fprintf('theta%d: %f \n',ii-1,rad2deg(thetas(ii)));
    end
    finalImageDistanceOurs = imageCoords(end,1)-opticalSystem(5,1);
    finalImageDistanceElagha = 17.768432;
    fprintf('Elegha gives a final image distance of %4.2f.\n',finalImageDistanceElagha)
    fprintf('We obtain: i5 - c5 = K4 = %4.2f \n',finalImageDistanceOurs);
    assert( abs(finalImageDistanceOurs - finalImageDistanceElagha) < 1e-4);
%}
%{
    %% Pupil through cornea
    % A model of the passage of a point on the pupil perimeter through
    % the axial cross-section of the cornea (units in mm)
    sceneGeometry = createSceneGeometry('aqueousRefractiveIndex',[]);
    outputRay = rayTraceCenteredSurfaces([sceneGeometry.eye.pupil.center(1) 2], deg2rad(-15), sceneGeometry.refraction.opticalSystem.p1p2, true);
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

% Pre-allocate our loop variables; set the values for at the first surface
% (initial position of ray)
aVals = zeros(nSurfaces);
aVals(1) = 1;
curvature = zeros(nSurfaces,1);
curvatureCenters = zeros(nSurfaces,1);
intersectionCoords=zeros(nSurfaces,2);
intersectionCoords(1,:)=coordsInitial;
thetas = zeros(nSurfaces,1);
thetas(1)=thetaInitial;
relativeIndices = zeros(nSurfaces,1);
relativeIndices(1) = 1;
imageCoords=zeros(nSurfaces,2);
% The initial image location is a projection of the initial ray back to the
% optical axis.
imageCoords(1,:)=[coordsInitial(1)-(coordsInitial(2)/tan(thetaInitial)) 0];

% Build the local optical system. Replace the center and radius of the
% first surface with the point of intersection of the initial ray with the
% optical axis, and set the radius to zero.
opticalSystem = zeros(nSurfaces,nDims);
opticalSystem(1,1)=imageCoords(1,1);
opticalSystem(1,2:size(opticalSystemIn,2)-1)=0;
opticalSystem(1,size(opticalSystemIn,2))=opticalSystemIn(1,end);
opticalSystem(2:end,:) = opticalSystemIn(2:end,:);
curvatureCenters(1) = opticalSystem(1,1);
% If a single radius value was passed for the surfaces, copy this value
% over to define the radius for these spheres along the horizontal axis
if size(opticalSystemIn,2)==3
    opticalSystem = [opticalSystem(:,1) opticalSystem(:,2) opticalSystem(:,2) opticalSystem(:,3)];
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
    % and the radius of curvature of the lens at that point
    [intersectionCoords(ii,:),curvature(ii)] = calcEllipseIntersect(intersectionCoords(ii-1,:), thetas(ii-1), opticalSystem(ii,1), opticalSystem(ii,2:3) );
    % Check if there if the ray was not intersecting
    if isnan(curvature(ii))
        warning('rayTraceCenteredSurfaces:nonIntersectingRay','Ray did not intersect surface %d. Returning.',ii);
        return
    end
    % Find the curvature center, which is the position along the optical
    % axis for this surface, given the curvature encountered by the ray.
    curvatureCenters(ii) = opticalSystem(ii,1)-opticalSystem(ii,2)+curvature(ii);
    % The distance between the center of curvature of the current lens
    % surface and the center of curvature of the prior lens surface
    d = curvatureCenters(ii)-curvatureCenters(ii-1);
    % The relative refractive index of the prior medium to the medium of
    % the surface that the ray is now impacting
    relativeIndices(ii)=opticalSystem(ii-1,end)/opticalSystem(ii,end);
    % Equation 54 of Elagha
    aVals(ii) = ...
        (1/curvature(ii))*(relativeIndices(ii-1).*aVals(ii-1).*curvature(ii-1)+d.*sin(thetas(ii-1)));
    % Check if the incidence angle is above the critical angle for the
    % relative refractive index at the surface interface.
    if abs((aVals(ii)*relativeIndices(ii))) > 1
        warning('rayTraceCenteredSurfaces:criticalAngle','Angle of incidence for surface %d greater than critical angle. Returning.',ii);
        return
    end
    % Find the angle of the ray after it enters the current surface
    thetas(ii) = thetas(ii-1) - asin(aVals(ii)) + asin(aVals(ii).*relativeIndices(ii));
    % Find the coordinates at which the ray, after making contact
    % with the current surface, would contact (or originate from)
    % the optical axis
    imageCoords(ii,:) = [-(intersectionCoords(ii,2)-tan(thetas(ii))*intersectionCoords(ii,1))/tan(thetas(ii)) 0];
    % Update the plot
    if figureFlag.show
        % add this lens surface
        if figureFlag.surfaces
            plotLensArc(opticalSystem(ii,:))
        end
        % plot the line for the virtual image
        if figureFlag.imageLines
            % Plot the prior virtual image line
            plot([imageCoords(ii-1,1) intersectionCoords(ii-1,1)],[imageCoords(ii-1,2) intersectionCoords(ii-1,2)],'--b');
        end
        % plot the line for the path of the ray
        if figureFlag.rayLines
            plot([intersectionCoords(ii-1,1) intersectionCoords(ii,1)],[intersectionCoords(ii-1,2) intersectionCoords(ii,2)],'-r');
        end
    end
end


%% Finish and clean up
% Assemble an output which is the unit vector for the final ray
slope = tan(thetas(nSurfaces)+pi);
norm = sqrt(slope^2+1);
outputRay = [intersectionCoords(nSurfaces,:); [intersectionCoords(nSurfaces,1)+(1/norm) intersectionCoords(nSurfaces,2)+(slope/norm)]];

% Complete the plot
if figureFlag.show
    % Plot the final virtual image
    if figureFlag.imageLines
        plot([imageCoords(nSurfaces,1) intersectionCoords(nSurfaces,1)],[imageCoords(nSurfaces,2) intersectionCoords(nSurfaces,2)],'-b');
    end
    % Plot the ray path, which is a triple-length output ray. Also, add a
    % mark to indicate the initial coordinates of the ray
    if figureFlag.rayLines
        slope = tan(thetas(nSurfaces)+pi);
        norm = sqrt(slope^2+1);
        finalRay = [intersectionCoords(nSurfaces,:); [intersectionCoords(nSurfaces,1)+(3/norm) intersectionCoords(nSurfaces,2)+(3*slope/norm)]];
        plot([finalRay(1,1) finalRay(2,1)],[finalRay(1,2) finalRay(2,2)],'-r');
        plot(intersectionCoords(1,1),intersectionCoords(1,2),'xr');
    end
    % Plot the output unit ray vector
    if figureFlag.finalUnitRay
        plot([outputRay(1,1) outputRay(2,1)],[outputRay(1,2) outputRay(2,2)],'-g');
    end
    % Replot the refline
    if figureFlag.refLine
        yl = [min([imageCoords(:,1); intersectionCoords(:,1); opticalSystem(:,1)]),max([imageCoords(:,1); intersectionCoords(:,1); opticalSystem(:,1)])];
        yl(1) = yl(1)-(yl(1)/3)*sign(yl(1));
        yl(2) = yl(2)+yl(1)/3*sign(yl(1));
        plot(yl,[0,0],'-k');
    end
    % Add some labels
    if figureFlag.textLabels
        yl = max(abs(intersectionCoords(:,2)));
        yPos = zeros(nSurfaces,1)-(yl/2);
        plot(opticalSystem(:,1),zeros(nSurfaces,1),'+k');
        labels = strseq('c',1:1:nSurfaces);
        text(opticalSystem(:,1),yPos,labels,'HorizontalAlignment','center');
        labels = strseq('i',1:1:nSurfaces);
        plot(imageCoords(:,1),zeros(nSurfaces,1),'*r');
        text(imageCoords(:,1),yPos,labels,'HorizontalAlignment','center');
    end
    % Add a legend
    if figureFlag.legend
        hSub = subplot(3,3,8);
        plot(nan, nan,'-r');
        hold on
        plot(nan, nan,'--b');
        plot(nan, nan,'-b');
        plot(nan, nan,'-g');
        set(hSub, 'Visible', 'off');
        legend({'ray path','virtual ray','final virtual ray','output unit ray vector'},'Location','north', 'Orientation','vertical');
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


function [intersectionCoords, curvature] = calcEllipseIntersect(coordsInitial, theta, ellipseCenterZ, ellipseRadii )
% Returns coords and curvature of an ellipse intersected by a ray
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
intersectionCoords = [nan, nan];
curvature = nan;

% Store the sign of the radius values. They radii must have the same sign.
if prod(sign(ellipseRadii))==-1
    error('rayTraceCenteredSurfaces:incompatibleConvexity','The radii of the elliptical lens surface must have the same sign.');
end
radiiSign = sign(ellipseRadii(1));
% Convert the radii to their absolute values
ellipseRadii = abs(ellipseRadii);

% Convert the ray position and theta to the slope (m) and y-axis
% intercept(c) of a line.
m = tan(theta);
c = coordsInitial(2)-m*coordsInitial(1);

% Obtain the pair of z and h coordinates at which the line will intersect
% the ellipse
M = (1/ellipseRadii(1)^2) + (1/ellipseRadii(2)^2)* m^2;
N = 2*(1/ellipseRadii(2)^2)*m*c - 2*(1/ellipseRadii(1)^2)*ellipseCenterZ;
O = (1/ellipseRadii(2)^2)*c^2 + (1/ellipseRadii(1)^2)*ellipseCenterZ^2 -1 ;
determinant = (N^2 - 4* M * O);

% Detect if we have a tangential or non-intersecting ray
if ~isreal(determinant) || determinant<0
    return
end

zCoords = [(- N + sqrt(determinant))/ (2*M), (- N - sqrt(determinant))/ (2*M)];
hCoords = [m*zCoords(1) + c, m*zCoords(2) + c];

% If the radiiSign is positive, report the coordinates on the left-hand
% side of the ellipse, otherwise report the coordinates on the right
if radiiSign<0
    [~,rightIdx] = max(zCoords);
    intersectionCoords = [zCoords(rightIdx),hCoords(rightIdx)];
else
    [~,leftIdx] = min(zCoords);
    intersectionCoords = [zCoords(leftIdx),hCoords(leftIdx)];
end

% Detect if we have a tangential or non-intersecting ray
P = (intersectionCoords(1)-ellipseCenterZ)/ellipseRadii(1);
if ~isreal(P)
    return
end
if P>1 || P<-1
    return
end

% Calculate the radius of curvature at the point of intersection. If
% radiiSign is negative, report a negative radius of curvature
t = acos(P);
curvature = radiiSign*((ellipseRadii(1)^2*sin(t)^2 + ellipseRadii(2)^2*cos(t)^2)^(3/2))/(ellipseRadii(1)*ellipseRadii(2));

end