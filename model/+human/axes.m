function axes = axes( eye )

%{
% Obtain the distance in mm between the optic disc and fovea.
[geoDistance,~,~,geodeticPathCoords] = quadric.panouGeodesicDistance(S,axes.visual.geodetic,axes.opticDisc.geodetic);

% Check if the returned distance value is nan
if isnan(geoDistance)
    % The geodesic crosses the umbilicus, and thus the Panou solution is
    % undefined.
    geodeticPathCoords = [];
    fprintf('Surface geodesic distance fovea -> optic disc center is undefined due to umbilical crossing\n');
else
    fprintf('Surface geodesic distance fovea -> optic disc center = %0.2f mm\n',geoDistance);
end

eucDistance = sqrt(sum((quadric.ellipsoidalGeoToCart(axes.visual.geodetic,S)-quadric.ellipsoidalGeoToCart(axes.opticDisc.geodetic,S)).^2));
fprintf('Euclidean distance fovea -> optic disc center = %0.2f mm\n',eucDistance);

% Plot the surface
figure
boundingBox = sceneGeometry.refraction.retinaToPupil.opticalSystem(1,12:17);
quadric.plotSurface(S,boundingBox,[0.9 0.9 0.9],0.8);
camlight
lighting gouraud
hold on

% Plot lines of constant beta
for beta = -90:10:50
    coords =[];
    for omega = -180:3:180
        coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
    end
    plot3(coords(:,1),coords(:,2),coords(:,3),'-b');
end

% Plot lines of constant omega
for omega = -180:10:180
    coords =[];
    for beta = -90:3:50
        coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
    end
    plot3(coords(:,1),coords(:,2),coords(:,3),'-g');
end

% Add the retinal landmarks
eyePoint = quadric.ellipsoidalGeoToCart( axes.optical.geodetic, S );
plot3(eyePoint(1),eyePoint(2),eyePoint(3),'+r','MarkerSize',10);

eyePointFovea = quadric.ellipsoidalGeoToCart( axes.visual.geodetic, S );
plot3(eyePointFovea(1),eyePointFovea(2),eyePointFovea(3),'*r','MarkerSize',10);

eyePointOpticDisc = quadric.ellipsoidalGeoToCart( opticDiscG, S );
plot3(eyePointOpticDisc(1),eyePointOpticDisc(2),eyePointOpticDisc(3),'*m','MarkerSize',10);

if ~isempty(geodeticPathCoords)
    plot3(geodeticPathCoords(:,1),geodeticPathCoords(:,2),geodeticPathCoords(:,3),'-.y');
else
    plot3([eyePointFovea(1) eyePointOpticDisc(1)],[eyePointFovea(2) eyePointOpticDisc(2)],[eyePointFovea(3) eyePointOpticDisc(3)],'-y');
end
%}


% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% Set some fmincon options we will be using below
options = optimoptions(@fmincon,...
    'Display','off');


%% optical axis
% Eye axes are specified as rotations (in degrees) within the eye
% world coordinate frame for azimuth, elevation, and rotation. Axes
% are defined relative to the optical axis, which itself is set to
% be aligned with the p1 dimension of the eye world coordinate
% frame.
axes.optical.degField = [0 0 0];
axes.optical.geodetic = [-90 -90 0];
axes.optical.cartesian = quadric.ellipsoidalGeoToCart(axes.optical.geodetic,S)';


%% visual axis
% Set the desired angle in degrees of visual field between the optic and
% visual axes of the eye (effectively, kappa). We assume an azimuth alpha
% of 5.8 degrees for an emmetropic right eye (Figure 8 of Mathur 2013). We
% assume an elevation alpha of 2.5 degrees.
switch eye.meta.eyeLaterality
    case 'Right'
        axes.visual.degField = [5.8 2.5 0];
    case 'Left'
        axes.visual.degField = [-5.8 2.5 0];
end

% We now calculate the location on the retina corresponding to this visual
% angle. The objective function is the difference in the visual field
% position of a candidate retinal point (G) and the desired value. Note
% that the elevational angle is inverted, to account for the reversed
% rotation convention that we have for positive rotation values elevating
% the eye.
myObj = @(G) sum(angdiff(deg2rad(visualAngleBetweenRetinalCoords(eye,axes.optical.geodetic,G)),deg2rad(axes.visual.degField(1:2).*[1 -1])).^2).*1e100;

% Supply an initial guess as to the geodetic coords of the fovea, guided by
% spherical ametropia and laterality.
if (-0.5<eye.meta.sphericalAmetropia) && (eye.meta.sphericalAmetropia<0.5)
    x0 = [-85 -121 0];
end
if isequal(eye.meta.eyeLaterality,'Left')
    x0(2)=-x0(2);
end

% Perform the search
axes.visual.geodetic = fmincon(myObj, x0, [], [], [], [], [-89 -180 0], [-70 0 0], [], options);

% Confirm that the visual angles are as desired
assert(max(abs(visualAngleBetweenRetinalCoords(eye,axes.optical.geodetic,axes.visual.geodetic)-axes.visual.degField(1:2).*[1 -1]))<1e-2);

% Obtain the Cartesian coordinates of the fovea
axes.visual.cartesian = quadric.ellipsoidalGeoToCart(axes.visual.geodetic,S)';


%% optic disc axis (physiologic blind spot)
% Values taken from Safren 1993 for their dim stimulus, under the
% assumption that this will be the most accurate given the minimization of
% light scatter.
switch eye.meta.eyeLaterality
    case 'Right'
        axes.opticDisc.degField = [-16.02 -1.84 0];
    case 'Left'
        axes.opticDisc.degField = [16.02 -1.84 0];
end

% Define the objective. Again note that the vertical target angle in
% degrees of visual field is reversed.
myObj = @(G) sum(angdiff(deg2rad(visualAngleBetweenRetinalCoords(eye,axes.visual.geodetic,G)),deg2rad(axes.opticDisc.degField(1:2).*[1 -1])).^2).*1e100;

% Define some x0 values to initialize the search, varying by ametropia
if isequal(quadric.dimensionSizeRank(S),[1 3 2])
    x0 = [-78 -5 0];
    lb = [-89 -30 0];
    ub = [-70 30 0];
end
if isequal(quadric.dimensionSizeRank(S),[1 2 3])
    x0 = [-75 120 0];
    lb = [-89 -30 0];
    ub = [-70 180 0];
end

% Perform the search
axes.opticDisc.geodetic = fmincon(myObj, x0, [], [], [], [], lb, ub, [], options);

% Confirm that the visual angles are as desired
assert(max(abs(visualAngleBetweenRetinalCoords(eye,axes.visual.geodetic,axes.opticDisc.geodetic)-axes.opticDisc.degField(1:2).*[1 -1]))<1e-1);

% Obtain the Cartesian coordinates of the optic disc
axes.opticDisc.cartesian = quadric.ellipsoidalGeoToCart(axes.opticDisc.geodetic,S)';


end

