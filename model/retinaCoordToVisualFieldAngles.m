function foo
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

eye = modelEyeParameters('sphericalAmetropia',0);
opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName', 'retinaToPupil' );
S = opticalSystem(1,1:10);
boundingBox = opticalSystem(1,12:17);



vitreousChamberApex = quadric.ellipsoidalGeoToCart( [-90 90 0], S );
eyePoint = quadric.ellipsoidalGeoToCart( [-85 85], S );




% Plot the surface
quadric.plotSurface(S,boundingBox,[0.9 0.9 0.9],0.8)
camlight
lighting gouraud
hold on
% Plot lines of constant beta
for beta = -90:10:90
    coords =[];
    for omega = -180:3:180
        coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
    end
    plot3(coords(:,1),coords(:,2),coords(:,3),'-b');
end
fprintf('Lines of constant beta in blue\n');
% Plot lines of constant omega
for omega = -180:10:180
    coords =[];
    for beta = -90:3:90
        coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
    end
    plot3(coords(:,1),coords(:,2),coords(:,3),'-g');
end
fprintf('Lines of constant omega in green\n');
plot3(eyePoint(1),eyePoint(2),eyePoint(3),'*y');
plot3(vitreousChamberApex(1),vitreousChamberApex(2),vitreousChamberApex(3),'*m');
plot3(geodeticPathCoords(:,1),geodeticPathCoords(:,2),geodeticPathCoords(:,3),'+r');


% args = {eye.pupil.center([2 3 1])', ...
%     eye.rotationCenters, ...
%     opticalSystem};

% eyePose = [0 0 0 2];
% [~, ~, angle_p1p2, angle_p1p3] = ...
%     virtualImageFunc(...
%     eyePoint', eyePose, args{:});


end %foo