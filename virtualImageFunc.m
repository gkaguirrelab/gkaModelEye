function [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, cameraPositionTranslation, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3 )
% Returns the virtual image location of a point in eyeWorld coordinates
%
% Syntax:
%  [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3 )
%
% Description:
%   This routine returns the virtual image location of a point that has
%   passed through an optical system. The function requires specification
%   of the eyeWorld coordinates of the point, the pose of the eye, as well
%   as several features of sceneGeometry. These last four input arguments
%   can be supplied with sceneGeometry.refraction.args{:}.
%
% Inputs:
%   eyeWorldPoint         - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The pupil
%                           radius value is unused by this routine.
%   cameraPositionTranslation - Equal to sceneGeometry.cameraPosition.
%                           translation
%   rotationCenters       - Equal to sceneGeometry.eye.rotationCenters
%   opticalSystem_p1p2, opticalSystem_p1p3 - Equal to sceneGeometry.
%                           refraction.opticalSystem.p1p2 and .p1p3.
%
% Outputs:
%   virtualEyeWorldPoint  - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   nodalPointIntersectError - The distance (in mm) between the nodal point
%                           of the camera within the Z camera plane and the
%                           point of intersection of the ray arising from
%                           the eyeWorld point.
%
% Examples:
%{
    % Basic example that finds the virtual image location for a point from
    % the top of a 2mm radius exit pupil, with the eye posed straight
    % ahead, and the camera in its default location.
    sceneGeometry = createSceneGeometry('forceMATLABVirtualImageFunc',true);
    % Assemble the args for the virtualImageFunc
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.opticalSystem.p1p2, ...
    	sceneGeometry.refraction.opticalSystem.p1p3};
    [virtualEyeWorldPoint, nodalPointIntersectError] = sceneGeometry.refraction.handle( [-3.7 2 0], [0 0 0 2], args{:} );
    % Test output against cached value
    virtualEyeWorldPointCached = [-3.700000000000000   2.260311180204983   0.000000000000000];
    assert(max(abs(virtualEyeWorldPoint - virtualEyeWorldPointCached)) < 1e-6)
%}
%{
    %% Confirm that nodalPointIntersectError remain small across eye poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Perform 100 forward projections with randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*45, (rand(nPoses,1)-0.5)*40, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    clear nodalPointIntersectError
    for pp = 1:nPoses
    	[~,~,~,~,~,nodalPointIntersectError(:,pp)]=pupilProjection_fwd(eyePoses(pp,:),sceneGeometry);
    end
    % Make sure the nodalPointIntersectError is small and not
    % systematically related to eyePose
    figure
    plot(sqrt(eyePoses(:,1).^2+eyePoses(:,2).^2),median(nodalPointIntersectError),'.r')
    xlabel('Euclidean rotation distance [deg]');
    ylabel('Ray trace nodal error [mm]');
%}

%% Handle ray trace warnings
coder.extrinsic('warning')
warnState = warning();
warning('off','rayTraceCenteredSurfaces:criticalAngle');
warning('off','rayTraceCenteredSurfaces:nonIntersectingRay');


%% Find the p1p2 and p1p3 thetas
% For this eyeWorld point, we find the theta value in the p1p2 and p1p3
% planes that results in a ray (after passing through the optical system)
% that minimizes the distance between the camera node and point of
% intersection with the camera plane

nodalPointIntersectError = Inf;
nodalErrorTolerance = 1e-4;
searchIterTolerance = 6;
theta_p1p2 = 1e-4;
theta_p1p3 = 1e-4;
searchingFlag = true;
searchIter = 0;

% Set fminsearch options to tolerate an error of 1e-2, and to make changes
% in theta as small as 1e-6.
% Perform the search
TolFun = 1e-2; % nodal point intersection error to tolerate
TolX = 1e-6; % precision with which theta is estimated
options = optimset('TolFun',TolFun,'TolX',TolX);

% Set a scoped variable that detects if we encountered a bad ray trace.
badTraceFlag = false;

while searchingFlag
    cameraNodeDistanceError3D_p1p2 = @(thetaX) calcCameraNodeDistanceError3D(eyeWorldPoint, thetaX, theta_p1p3, eyePose, cameraPositionTranslation, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3);
    [theta_p1p2,~]=fminsearch(@myObj_p1p2,theta_p1p2,options);
    % Detect if we hit a bad ray trace. We will tolerate a bad trace if the
    % solution theta is very close to zero (within an order of magnitude of
    % TolX). This is because the ray tracing solution is undefined exactly at
    % zero, so bad ray traces occur when we search for thetas close to zero. If
    % we have encountered a bad ray trace away from zero, exit the routine.
    if badTraceFlag && theta_p1p2 > (TolX*10)
        virtualEyeWorldPoint = nan(1,3);
        nodalPointIntersectError = Inf;
        % Restore the warning state
        warning(warnState);
        return
    else
        badTraceFlag = false;
    end
    cameraNodeDistanceError3D_p1p3 = @(thetaX) calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, thetaX, eyePose, cameraPositionTranslation, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3);
    [theta_p1p3,nodalPointIntersectError]=fminsearch(@myObj_p1p3,theta_p1p3,options);
    % Detect if we hit a bad ray trace.
    if badTraceFlag && theta_p1p3 > (TolX*10)
        virtualEyeWorldPoint = nan(1,3);
        nodalPointIntersectError = Inf;
        % Restore the warning state
        warning(warnState);
        return
    else
        badTraceFlag = false;
    end
    % Determine if we have reached one of our stopping criteria
    searchIter = searchIter+1;
    if nodalPointIntersectError<nodalErrorTolerance || searchIter>searchIterTolerance
        searchingFlag = false;
    end
end
    function fval = myObj_p1p2(x)
        fval = cameraNodeDistanceError3D_p1p2(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end
    function fval = myObj_p1p3(x)
        fval = cameraNodeDistanceError3D_p1p3(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end


%% Obtain the virtual image location
% With both theta values calculated, now obtain the virtual image
% ray arising from the pupil plane that reflects the corneal optics
virtualImageRay = calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, opticalSystem_p1p2, opticalSystem_p1p3);

% Restore the warning state
warning(warnState);

% Extract the origin of the ray, which is the virtual image eyeWorld point
virtualEyeWorldPoint = virtualImageRay(1,:);

end % virtualImageFunc -- MAIN



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% calcCameraNodeDistanceError3D
function distance = calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, cameraPositionTranslation, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3)
% 3D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3)
%
% Description:
%   This function is similar to the calcCameraNodeDistanceError2D_p1p2,
%   except that it takes as input theta in both the p1p2 and p1p3 planes.
%   The distance value that is returned is still the Euclidean distance
%   between the intersection point of the output ray on the Z camera plane
%   and the nodal point of the camera.


% Ray trace for these thetas
outputRayEyeWorld2D_p1p2 = rayTraceCenteredSurfaces([eyeWorldPoint(1), eyeWorldPoint(2)], theta_p1p2, opticalSystem_p1p2);
outputRayEyeWorld2D_p1p3 = rayTraceCenteredSurfaces([eyeWorldPoint(1), eyeWorldPoint(3)], theta_p1p3, opticalSystem_p1p3);

% If we received a ray-trace error, then return Inf for the distance
if isempty(outputRayEyeWorld2D_p1p2) || isempty(outputRayEyeWorld2D_p1p3)
    distance = Inf;
    return
end

% Shift the p1p3 ray to have the same initial p1 value as the p1p2 ray
slope =(outputRayEyeWorld2D_p1p2(2,2)-outputRayEyeWorld2D_p1p2(1,2))/(outputRayEyeWorld2D_p1p2(2,1)-outputRayEyeWorld2D_p1p2(1,1));
zOffset=outputRayEyeWorld2D_p1p2(1,1)-eyeWorldPoint(1);
outputRayEyeWorld2D_p1p2(:,1)=outputRayEyeWorld2D_p1p2(:,1)-zOffset;
outputRayEyeWorld2D_p1p2(:,2)=outputRayEyeWorld2D_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld2D_p1p3(2,2)-outputRayEyeWorld2D_p1p3(1,2))/(outputRayEyeWorld2D_p1p3(2,1)-outputRayEyeWorld2D_p1p3(1,1));
zOffset=outputRayEyeWorld2D_p1p3(1,1)-eyeWorldPoint(1);
outputRayEyeWorld2D_p1p3(:,1)=outputRayEyeWorld2D_p1p3(:,1)-zOffset;
outputRayEyeWorld2D_p1p3(:,2)=outputRayEyeWorld2D_p1p3(:,2)-(zOffset*slope);

% Combine into a single, 3D ray
outputRayEyeWorld3D=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) outputRayEyeWorld2D_p1p3(1,2);...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) outputRayEyeWorld2D_p1p3(2,2)];

% prepare to rotate the outputRay into the sceneWorld coordinates
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(-eyePose(2)) 0 sind(-eyePose(2)); 0 1 0; -sind(-eyePose(2)) 0 cosd(-eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Copy over the outputRay from eye to head world
outputRayHeadWorld3D=outputRayEyeWorld3D;

% Torsion
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenters.tor(dim);
    end
end
outputRayHeadWorld3D = (RotTor*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenters.tor(dim);
    end
end
% Elevation
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenters.ele(dim);
    end
end
outputRayHeadWorld3D = (RotEle*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenters.ele(dim);
    end
end
% Azimuth
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenters.azi(dim);
    end
end
outputRayHeadWorld3D = (RotAzi*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenters.azi(dim);
    end
end

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld3D = outputRayHeadWorld3D(:,[2 3 1]);

% Obtain an expression for X and Y distances between the nodal point of the
% camera in the sceneWorld plane and the point at which the ray will strike
% the plane that contains the camera
slope_xZ =(outputRaySceneWorld3D(2,1)-outputRaySceneWorld3D(1,1))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
slope_yZ =(outputRaySceneWorld3D(2,2)-outputRaySceneWorld3D(1,2))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
cameraPlaneX = outputRaySceneWorld3D(1,1)+((cameraPositionTranslation(3)-outputRaySceneWorld3D(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld3D(1,2)+((cameraPositionTranslation(3)-outputRaySceneWorld3D(1,3))*slope_yZ);

% Compute the Euclidean distance between the point of intersection and the
% nodal point of the camera.
distance = sqrt((cameraPositionTranslation(1)-cameraPlaneX)^2 + ...
        (cameraPositionTranslation(2)-cameraPlaneY)^2 );

end % calcCameraNodeDistanceError3D



%% calcVirtualImageRay
function [outputRayEyeWorld3D] = calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, opticalSystem_p1p2, opticalSystem_p1p3)
% Returns the unit vector virtual image ray for the initial depth position
%
% Syntax:
%  [outputRayEyeWorld] = calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, opticalSystem_p1p2, opticalSystem_p1p3)
%
% Description:
%   For a given point in eyeWorld coordinates, and for a given pair of
%   theta values, this function returns the ray that corresponds to the
%   virtual image that arises from the optical system, with the initial
%   point of the virtual image being at the same p1 position as the
%   object point.
%
%   Practically, once the p1p2 and p1p3 thetas are found, this function is
%   used to obtain the position within the eyeWorld coordinate frame that
%   is the apparent location of the point after refraction through the
%   cornea. For this reason, only first row of values are used by the
%   calling function.
%
% Inputs:
%   eyeWorldPoint
%   theta_p1p2, theta_p1p3 - Scalar in units of radians. The angle WRT the
%                           optical axis of the initial ray. The function
%                           is undefined for theta = 0 (i.e., a paraxial
%                           ray) and will return nan. Also, absolute values
%                           of pi correspond to a vertical ray that would
%                           not intersect with the optical system and thus
%                           will return nan. There are other combinations
%                           of eyeWorld positions and thetas that will
%                           return nan given the particular path of the ray
%                           through the optical system.
%
% Outputs:
%   outputRay             - A 2x3 matrix that is the unit vector of a ray
%                           in eyeWorld coordinates. outputRay(1,:) is
%                           the origin point of the ray, corresponding
%                           to dimensions [p1 p2 p3], and the value of
%                           p1 set to be equal to the input value of p1.
%                           The values in outputRay(2,:) give the position
%                           of the unit vector.
%


% Ray trace for these thetas
outputRayEyeWorld2D_p1p2 = rayTraceCenteredSurfaces([eyeWorldPoint(1), eyeWorldPoint(2)], theta_p1p2, opticalSystem_p1p2);
outputRayEyeWorld2D_p1p3 = rayTraceCenteredSurfaces([eyeWorldPoint(1), eyeWorldPoint(3)], theta_p1p3, opticalSystem_p1p3);

% If we received a ray-trace error, then return nans for output ray
if isempty(outputRayEyeWorld2D_p1p2) || isempty(outputRayEyeWorld2D_p1p3)
    outputRayEyeWorld3D = nan(2,3);
    return
end

% Adjust the p1 (optical axis) position of the rays to have their initial
% position at the same p1
slope =(outputRayEyeWorld2D_p1p2(2,2)-outputRayEyeWorld2D_p1p2(1,2))/(outputRayEyeWorld2D_p1p2(2,1)-outputRayEyeWorld2D_p1p2(1,1));
zOffset=outputRayEyeWorld2D_p1p2(1,1)-eyeWorldPoint(1);
outputRayEyeWorld2D_p1p2(:,1)=outputRayEyeWorld2D_p1p2(:,1)-zOffset;
outputRayEyeWorld2D_p1p2(:,2)=outputRayEyeWorld2D_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld2D_p1p3(2,2)-outputRayEyeWorld2D_p1p3(1,2))/(outputRayEyeWorld2D_p1p3(2,1)-outputRayEyeWorld2D_p1p3(1,1));
zOffset=outputRayEyeWorld2D_p1p3(1,1)-eyeWorldPoint(1);
outputRayEyeWorld2D_p1p3(:,1)=outputRayEyeWorld2D_p1p3(:,1)-zOffset;
outputRayEyeWorld2D_p1p3(:,2)=outputRayEyeWorld2D_p1p3(:,2)-(zOffset*slope);

outputRayEyeWorld3D = zeros(2,3);

% Combine the two dimensions into a single, 3D ray
outputRayEyeWorld3D(1,:) = [outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) outputRayEyeWorld2D_p1p3(1,2)];
outputRayEyeWorld3D(2,:) = [outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) outputRayEyeWorld2D_p1p3(2,2)];

end % calcVirtualImageRay


