function [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3 )
% Returns the virtual image location of a point in eyeWorld coordinates
%
% Syntax:
%  [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3 )
%
% Description:
%   This routine returns the virtual image location of a point that has passed through an
%   optical system. The computation is assembled step-wise:
%       traceOpticalSystemFuncHandle - 2D ray tracing through the cornea and any
%           corrective lenses
%       calcCameraNodeDistanceError2D - 2D distance of ray intersection on 
%           camera plane from camera node
%       calcCameraNodeDistanceError3D - 3D distance of ray intersection on 
%           camera plane from camera node
%       calcVirtualImageRay - Returns the unit vector virtual image ray for
%           the initial depth position
%
% Inputs:
%   eyeWorldPoint         - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius in mm.
%   opticalSystem         - Equal to sceneGeometry.opticalSystem
%   extrinsicTranslationVector - Equal to sceneGeometry.
%                               extrinsicTranslationVector
%   rotationCenters       - Equal to sceneGeometry.eye.rotationCenters
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
    sceneGeometry = createSceneGeometry();
    [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( [-3.7 2 0], [0 0 0 2], sceneGeometry.virtualImageFunc.args{:} );
    % Test output against value computed on April 10, 2018
    virtualEyeWorldPointStored = [-3.7000    2.2553    0.0000];
    assert(max(abs(virtualEyeWorldPoint - virtualEyeWorldPointStored)) < 1e-4)
%}


%% Handle ray trace warnings
coder.extrinsic('warning')
warnState = warning();
warning('off','rayTraceCenteredSurfaces:criticalAngle');
warning('off','rayTraceCenteredSurfaces:nonIntersectingRay');


%% Find the p1p2 theta
% For this eyeWorld point, we find the theta value in the p1p2 plane that
% results in a ray (after passing through the optical system) that
% minimizes the distance between the camera node and point of intersection
% with the camera plane

% Define an error function that is the distance between the nodal point of
% the camera and the point at which a ray intersects the plane that
% contains the camera, with the ray departing from the eyeWorld point at
% angle theta in the p1p2 plane.
cameraNodeDistanceError2D_p1p2 = @(theta_p1p2) calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2);

% Conduct an fminsearch to find the p1p2 theta that minimizes the node
% distance error. The options are adjusted to tolerate an error of 1e-2,
% and to make changes in theta as small as 1e-6. Because the errorFunc
% returns nan for values very close to zero, we initialize the search with
% a slightly non-zero value (1e-4)

% Set a scoped variable that detects if we encountered a bad ray trace.
badTraceFlag = false;

% Perform the search
options = optimset('TolFun',1e-2,'TolX',1e-6);
theta_p1p2=fminsearch(@myObj_p1p2,1e-4,options);
    function fval = myObj_p1p2(x)
        fval = cameraNodeDistanceError2D_p1p2(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end

% If we hit a bad ray trace, exit the routine
if badTraceFlag
    virtualEyeWorldPoint = nan(1,3);
    nodalPointIntersectError = Inf;
    % Restore the warning state
    warning(warnState);
    return
end

%% Find the p1p3 theta
% Given this p1p2 theta, we now find the p1p3 theta that further reduces
% the distance to the nodal point.
cameraNodeDistanceError3D = @(theta_p1p3) calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3);

% The fVal at the solution is the the total error (in mm) in both
% dimensions for intersecting the nodal point of the camera.

% Set a scoped variable that detects if we encountered a bad ray trace.
badTraceFlag = false;

% Perform the search
options = optimset('TolFun',1e-2,'TolX',1e-6);
[theta_p1p3, nodalPointIntersectError]=fminsearch(@myObj_3D,1e-4,options);
    function fval = myObj_3D(x)
        fval = cameraNodeDistanceError3D(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end

% If we hit a bad ray trace, exit the routine
if badTraceFlag
    virtualEyeWorldPoint = nan(1,3);
    nodalPointIntersectError = Inf;
    % Restore the warning state
    warning(warnState);
    return
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



%% calcCameraNodeDistanceError2D_p1p2
function distance = calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2)
% 2D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2)
%
% Description:
%   This function returns the distance between the nodal point of the
%   camera and a ray that has exited from the optical system. This distance
%   is calculated within the sceneWorld coordinates on an X-Y plane that is
%   positioned at the Z location of the camera. The point of intersection
%   of the ray upon the plane is found, and then the Euclidean distance
%   between this impact point and the nodal point of the camera is
%   returned.
%
%   This function is used to find a theta in the p1p2 plane that minimize
%   the distance between the the intersection point of the ray in the
%   camera plane and the nodal point of the camera. At a distance of zero,
%   the ray would enter the pin hole aperture of the camera and thus
%   produce a point on the resulting image. A subsequent search,
%   constrained by the p1p2 results, finds the p1p3 theta value.
%
% Inputs:
%   eyeWorldPoint
%   theta_p1p2            - Scalar in radians. The angle w.r.t. the optical
%                           axis of the initial ray. The function is
%                           undefined for theta = 0 (i.e., a paraxial ray)
%                           and will return nan. Also, absolute values of
%                           pi correspond to a vertical ray that would not
%                           intersect with the optical system and thus will
%                           return nan. There are other combinations of
%                           eyeWorld positions and thetas that will return
%                           nan given the particular path of the ray
%                           through the optical system.
%   eyePose
%   rotationCenters
%   traceOpticalSystemFuncHandle
%
% Outputs:
%   distance              - Scalar in units of mm. The Euclidean distance
%                           of the intersection point of the ray on the
%                           Z camera plane from the nodal point of the
%                           camera.
%


% Ray trace for this theta
outputRayEyeWorld2D_p1p2 = rayTraceCenteredSurfaces([eyeWorldPoint(1),eyeWorldPoint(2)],theta_p1p2, opticalSystem_p1p2);

% If we received a ray-trace error, then return Inf for the distance
if isempty(outputRayEyeWorld2D_p1p2)
    distance = Inf;
    return
end

% Add the p3 dimension
outputRayEyeWorld_p1p2=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) eyeWorldPoint(3);...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) eyeWorldPoint(3)];

% Prepare to rotate the outputRay into the sceneWorld coordinates
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(eyePose(2)) 0 sind(eyePose(2)); 0 1 0; -sind(eyePose(2)) 0 cosd(eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Copy eyeWorld rays over the HeadWorld variables
outputRayHeadWorld_p1p2 = outputRayEyeWorld_p1p2;

% For each of the two coordinates in each ray, shift the eyeWorld ray to
% the rotational center of the eye, rotate for this eye pose, then undo the
% centering

% Torsion
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenters.tor(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotTor*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenters.tor(dim);
    end
end
% Elevation
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenters.ele(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotEle*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenters.ele(dim);
    end
end
% Azimuth
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenters.azi(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotAzi*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenters.azi(dim);
    end
end

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld_p1p2 = outputRayHeadWorld_p1p2(:,[2 3 1]);

% Obtain an expression for X and Y distances between the nodal point of the
% camera in the sceneWorld plane and the point at which the ray will strike
% the plane that contains the camera
slope_xZ =(outputRaySceneWorld_p1p2(2,1)-outputRaySceneWorld_p1p2(1,1))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
slope_yZ =(outputRaySceneWorld_p1p2(2,2)-outputRaySceneWorld_p1p2(1,2))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
cameraPlaneX = outputRaySceneWorld_p1p2(1,1)+((extrinsicTranslationVector(3)-outputRaySceneWorld_p1p2(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p2(1,2)+((extrinsicTranslationVector(3)-outputRaySceneWorld_p1p2(1,3))*slope_yZ);

% Compute the Euclidean distance between the point of intersection and the
% nodal point of the camera.
distance = sqrt((extrinsicTranslationVector(1)-cameraPlaneX)^2 + ...
        (extrinsicTranslationVector(2)-cameraPlaneY)^2 );

end % calcCameraNodeDistanceError2D_p1p2



%% calcCameraNodeDistanceError3D
function distance = calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, extrinsicTranslationVector, rotationCenters, opticalSystem_p1p2, opticalSystem_p1p3)
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
RotEle = [cosd(eyePose(2)) 0 sind(eyePose(2)); 0 1 0; -sind(eyePose(2)) 0 cosd(eyePose(2))];
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
cameraPlaneX = outputRaySceneWorld3D(1,1)+((extrinsicTranslationVector(3)-outputRaySceneWorld3D(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld3D(1,2)+((extrinsicTranslationVector(3)-outputRaySceneWorld3D(1,3))*slope_yZ);

% Compute the Euclidean distance between the point of intersection and the
% nodal point of the camera.
distance = sqrt((extrinsicTranslationVector(1)-cameraPlaneX)^2 + ...
        (extrinsicTranslationVector(2)-cameraPlaneY)^2 );

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



