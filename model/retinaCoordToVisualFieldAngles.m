function foo
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sceneGeometry = createSceneGeometry('sphericalAmetropia',-1);

% Obtain the quadric form of the retinal surface
S = sceneGeometry.refraction.retinaToPupil.opticalSystem(1,1:10);

% Set some fmincon options we will be using below
options = optimoptions(@fmincon,...
    'Display','off');

% Define the location of the retinal apex in ellipsoidal geodetic coords.
retinalApexG = [-90; -90; 0];

% Set the desired angle in degrees of visual field between the optic and
% visual axes of the eye (effectively, kappa). Define the objective.
targetAngles = [5.8 -2.5];
myObj = @(G) sum(angdiff(deg2rad(visualAngleBetweenRetinalCoords(sceneGeometry,retinalApexG,G)),deg2rad(targetAngles)).^2).*1e100;

% Supply an initial guess as to the geodetic coords of the fovea
x0 = [-85;-110;0];

% Perform the search
foveaG = fmincon(myObj, x0, [], [], [], [], [-89;-180;0], [-70;0;0], [], options);

% Confirm that the visual angles are as desired
assert(max(abs(visualAngleBetweenRetinalCoords(sceneGeometry,retinalApexG,foveaG)-targetAngles))<1e-2);

% Set the desired angle in degrees of visual field between the visual axis
% of the eye and the center of the blind spot
targetAngles = [-16.02 1.84];
myObj = @(G) sum(angdiff(deg2rad(visualAngleBetweenRetinalCoords(sceneGeometry,foveaG,G)),deg2rad(targetAngles)).^2).*1e100;

% Start the search at the retinal apex.
if isequal(quadric.dimensionSizeRank(S),[1 3 2])
    x0 = [-78;-5;0];
    lb = [-89;-30;0];
    ub = [-70;30;0];
end
if isequal(quadric.dimensionSizeRank(S),[1 2 3])
    x0 = [-75;120;0];
    lb = [-89;-30;0];
    ub = [-70;180;0];
end

% Perform the search
opticDiscG = fmincon(myObj, x0, [], [], [], [], lb, ub, [], options);

% Confirm that the visual angles are as desired
assert(max(abs(visualAngleBetweenRetinalCoords(sceneGeometry,foveaG,opticDiscG)-targetAngles))<1e-1);

% Obtain the distance in mm between the optic disc and fovea.
[geoDistance,~,~,geodeticPathCoords] = quadric.panouGeodesicDistance(S,foveaG,opticDiscG);

% Check if the returned distance value is nan
if isnan(geoDistance)
    % The geodesic crosses the umbilicus, and thus the Panou solution is
    % undefined. In this case, I approximate the solution by taking the
    % geodesic on a bi-axial ellipsoid.
    radii = quadric.radii(quadric.alignAxes(S));
    a = mean(radii(3)); % The semi-major axis
    b = radii(1); % The semi-minor axis
    sharedLongitude=mean(abs([foveaG(2) opticDiscG(2)]));    
    geoDistanceA = distance('gc',foveaG(1),foveaG(2),-90,-sharedLongitude,[a sqrt(a^2-b^2)/a],'degrees');
    geoDistanceB = distance('gc',-90,sharedLongitude,opticDiscG(1),opticDiscG(2),[a sqrt(a^2-b^2)/a],'degrees');
    geoDistance = geoDistanceA+geoDistanceB;
    geodeticPathCoords = [];
    [lattrk,lontrk] = track('gc',[foveaG(1) -90],[foveaG(2) -sharedLongitude],[a sqrt(a^2-b^2)/a],'degrees',10);
    for ii=1:10
        geodeticPathCoords(end+1,:)=quadric.ellipsoidalGeoToCart([lattrk(ii);lontrk(ii);0],S);
    end
    [lattrk,lontrk] = track('gc',[opticDiscG(1) -90],[opticDiscG(2) sharedLongitude ],[a sqrt(a^2-b^2)/a],'degrees',10);
    for ii=1:10
        geodeticPathCoords(end+1,:)=quadric.ellipsoidalGeoToCart([lattrk(ii);lontrk(ii);0],S);
    end
    fprintf('Bi-axial estimated distance fovea -> optic disc center = %0.2f mm\n',geoDistance);
else
    fprintf('Distance fovea -> optic disc center = %0.2f mm\n',geoDistance);
end

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
eyePoint = quadric.ellipsoidalGeoToCart( retinalApexG, S );
plot3(eyePoint(1),eyePoint(2),eyePoint(3),'+r','MarkerSize',10);

eyePoint = quadric.ellipsoidalGeoToCart( foveaG, S );
plot3(eyePoint(1),eyePoint(2),eyePoint(3),'*r','MarkerSize',10);

eyePoint = quadric.ellipsoidalGeoToCart( opticDiscG, S );
plot3(eyePoint(1),eyePoint(2),eyePoint(3),'*m','MarkerSize',10);

plot3(geodeticPathCoords(:,1),geodeticPathCoords(:,2),geodeticPathCoords(:,3),'-y');


end %foo
