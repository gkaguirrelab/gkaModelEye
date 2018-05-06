function [virtualEyePoint, nodalPointIntersectError] = virtualImageFunc( eyePoint, eyePose, cameraTranslation, rotationCenters, opticalSystem )
% Returns the virtual image location of a point in eyeWorld coordinates
%
% Syntax:
%  [virtualEyePoint, nodalPointIntersectError] = virtualImageFunc( eyePoint, eyePose, cameraTranslation, rotationCenters, opticalSystem )
%
% Description:
%   This routine returns the virtual image location of a point that has
%   passed through an optical system. The function requires specification
%   of the eyeWorld coordinates of the point, the pose of the eye, as well
%   as several features of sceneGeometry.
%
% Inputs:
%   eyePoint              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The pupil
%                           radius value is unused by this routine.
%   cameraTranslation     - Equal to sceneGeometry.cameraPosition.
%                           translation
%   rotationCenters       - Equal to sceneGeometry.eye.rotationCenters
%   opticalSystem         - Equal to sceneGeometry.refraction.opticalSystem
%
% Outputs:
%   virtualEyeWorldPoint  - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   nodalPointIntersectError - The distance (in mm) between the nodal point
%                           of the camera and ray arising from the eyeWorld
%                           point after it exited the optical system.
%
% Examples:
%{
    % Basic example that finds the virtual image location for a point from
    % the top of a 2 mm radius exit pupil, with the eye posed straight
    % ahead, and the camera in its default location.
    sceneGeometry = createSceneGeometry('forceMATLABVirtualImageFunc',true);
    % Assemble the args for the virtualImageFunc
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.opticalSystem};
    [virtualEyePoint, nodalPointIntersectError] = sceneGeometry.refraction.handle( [sceneGeometry.eye.pupil.center(1) 2 0], [0 0 0 2], args{:} );
    % Test output against cached value
    virtualEyePointCached = [-4.000000000000000   2.290220465113814   0.000000000000000];
    assert(max(abs(virtualEyePoint - virtualEyePointCached)) < 1e-6)
%}
%{
    %% Confirm that nodalPointIntersectError remains small across eye poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Perform 100 forward projections with randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*60, (rand(nPoses,1)-0.5)*60, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
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


%% Find the p1p2 and p1p3 angles
% For this eyeWorld point, we find the angles of origin of a ray in the
% p1p2 and p1p3 planes that results in a ray (after passing through the
% optical system) that passes as close as possible to the camera node

% Pre-define the output variables to keep the compiler happy
nodalPointIntersectError = Inf;
virtualEyePoint = nan(1,3);

% Set some parameters for the search
nodalErrorTolerance = 1e-5;
searchIterTolerance = 6;
searchIter = 0;
searchingFlag = true;
angle_p1p2 = 1e-4;
angle_p1p3 = 1e-4;

% Set fminsearch options to tolerate an error of 1e-2, and to make changes
% in theta as small as 1e-6.
TolFun = 1e-2; % nodal point intersection error to tolerate
TolX = 1e-6; % precision with which theta is estimated
options = optimset('TolFun',TolFun,'TolX',TolX);

% Set a scoped variable that detects if we encountered a bad ray trace.
badTraceFlag = false;

% Enter a while loop that iteratively refines the theta values until
% criteria are met.
while searchingFlag

    % The distance error function for searching across p1p2 theta values
    cameraNodeDistanceError_p1p2 = @(angleX) calcCameraNodeDistanceError(eyePoint, angleX, angle_p1p3, eyePose, cameraTranslation, rotationCenters, opticalSystem);

    % Set the x0 value for the search, handling the issue that an initial
    % value of exactly zero breaks the search.
    if abs(angle_p1p2) < 1e-4
        x0 = 1e-4;
    else
        x0 = angle_p1p2;
    end
    % Perform the search
    [angle_p1p2,~]=fminsearch(@myObj_p1p2,x0,options);

    % Detect if we hit a bad ray trace. We will tolerate a bad trace if the
    % solution theta is very close to zero (within an order of magnitude of
    % TolX). This is because the ray tracing solution is undefined exactly
    % at zero, so bad ray traces occur when we search for thetas close to
    % zero. If we have encountered a bad ray trace away from zero, exit the
    % routine.
    if badTraceFlag && angle_p1p2 > (TolX*10)
        virtualEyePoint = nan(1,3);
        nodalPointIntersectError = Inf;
        % Restore the warning state
        warning(warnState);
        return
    else
        badTraceFlag = false;
    end
    
    % The distance error function for searching across p1p3 theta values    
    cameraNodeDistanceError_p1p3 = @(angleX) calcCameraNodeDistanceError(eyePoint, angle_p1p2, angleX, eyePose, cameraTranslation, rotationCenters, opticalSystem);

    % Set the x0 value for the search, handling the issue that an initial
    % value of exactly zero breaks the search.
    if abs(angle_p1p3) < 1e-4
        x0 = 1e-4;
    else
        x0 = angle_p1p3;
    end
    % Perform the search
    [angle_p1p3,nodalPointIntersectError]=fminsearch(@myObj_p1p3,x0,options);

    % Detect if we hit a bad ray trace.
    if badTraceFlag && angle_p1p3 > (TolX*10)
        virtualEyePoint = nan(1,3);
        nodalPointIntersectError = Inf;
        % Restore the warning state
        warning(warnState);
        return
    else
        badTraceFlag = false;
    end
    
    % Iterate the search count
    searchIter = searchIter+1;

    % Determine if we have met a stopping criterion
    if nodalPointIntersectError<nodalErrorTolerance || searchIter>searchIterTolerance
        searchingFlag = false;
    end
end
    % Nested objective functions. Needed to handle the badTrace behavior.
    function fval = myObj_p1p2(x)
        fval = cameraNodeDistanceError_p1p2(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end
    function fval = myObj_p1p3(x)
        fval = cameraNodeDistanceError_p1p3(x);
        if isinf(fval)
            badTraceFlag = true;
        end
    end


%% Obtain the virtual image location
% With both theta values calculated, now obtain the virtual image
% ray arising from the pupil plane that reflects the corneal optics
virtualEyePoint = calcVirtualImagePoint(eyePoint, angle_p1p2, angle_p1p3, opticalSystem);

% Restore the warning state
warning(warnState);

end % virtualImageFunc -- MAIN




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calcCameraNodeDistanceError
function distance = calcCameraNodeDistanceError(eyePoint, angle_p1p2, angle_p1p3, eyePose, cameraTranslation, rotationCenters, opticalSystem)
% Distance of ray intersection point on a camera plane from camera node
%
% Syntax:
%  distance = calcCameraNodeDistanceError(eyePoint, angle_p1p2, angle_p1p3, eyePose, cameraTranslation, rotationCenters, opticalSystem)
%
% Description:
%   This function returns the Euclidean distance between the nodal point of
%   the camera and a ray that has exited from the optical system of a
%   rotated eye. This distance is calculated within the eye coordinates.
%
%   This function is used to find angles in the p1p2 and p1p3 planes that
%   minimize the distance between the intersection point of the ray in the
%   camera plane and the nodal point of the camera. At a distance of zero,
%   the ray would enter the pin hole aperture of the camera and thus
%   produce a point on the resulting image.
%
% Inputs:
%   eyePoint
%   angle_p1p2, angle_p1p3 - Scalar in radians. The angle w.r.t. the optical
%                           axis of the initial ray. 
%   eyePose
%   cameraTranslation
%   rotationCenters
%   opticalSystem
%
% Outputs:
%   distance              - Scalar in units of mm. The Euclidean distance
%                           of the minimum distance between the nodal point
%                           of the camera and a ray exiting the optical
%                           system of the rotated eye. Set to Inf if an
%                           error is returned by rayTraceEllipsoids.
%


% Ray trace
outputRayEyeWorld = rayTraceEllipsoids(eyePoint, [angle_p1p2, angle_p1p3], opticalSystem);

% If we received a ray trace error, then return Inf for the distance
if isempty(outputRayEyeWorld)
    distance = Inf;
    return
end

% Counter-rotate the camera position point (in eyeWorld coordinate space)
% so that it is in a position w.r.t. the eye that is equivalent to if the
% eye had rotated
cameraRot = -eyePose;
RotAzi = [cosd(cameraRot(1)) -sind(cameraRot(1)) 0; sind(cameraRot(1)) cosd(cameraRot(1)) 0; 0 0 1];
RotEle = [cosd(-cameraRot(2)) 0 sind(-cameraRot(2)); 0 1 0; -sind(-cameraRot(2)) 0 cosd(-cameraRot(2))];
RotTor = [1 0 0; 0 cosd(cameraRot(3)) -sind(cameraRot(3)); 0 sind(cameraRot(3)) cosd(cameraRot(3))];

% Rearrange the camera translation dimensions to switch from world to eye
% coordinate space.
Pc = cameraTranslation([3 1 2])';

% Torsion
Pc=Pc-rotationCenters.tor;
Pc = (RotTor*Pc')';
Pc=Pc+rotationCenters.tor;
% Elevation
Pc=Pc-rotationCenters.ele;
Pc = (RotEle*Pc')';
Pc=Pc+rotationCenters.ele;
% Azimuth
Pc=Pc-rotationCenters.azi;
Pc = (RotAzi*Pc')';
Pc=Pc+rotationCenters.azi;

% Calculate the distance between the closest approach of the outputRay to
% the camera nodal point.
d = distancePointLine3d(Pc, [outputRayEyeWorld(1,:) (outputRayEyeWorld(2,:)-outputRayEyeWorld(1,:))]);

% Obtain the Euclidean distance in the 3 dimensions.
distance = sqrt(sum(d.^2));

end % calcCameraNodeDistanceError



%% calcVirtualImagePoint
function [virtualEyePoint] = calcVirtualImagePoint(eyePoint, angle_p1p2, angle_p1p3, opticalSystem)
% Returns the unit vector virtual image ray for the initial depth position
%
% Syntax:
%  [virtualEyePoint] = calcVirtualImagePoint(eyePoint, angle_p1p2, angle_p1p3, opticalSystem)
%
% Description:
%   For a given point in eyeWorld coordinates, and for a given pair of
%   angle values, this function returns the ray that corresponds to the
%   virtual image that arises from the optical system, with the initial
%   point of the virtual image being at the same p1 position as the
%   object point.
%
%   Practically, once the p1p2 and p1p3 angles are found, this function is
%   used to obtain the position within the eyeWorld coordinate frame that
%   is the apparent location of the point after refraction through the
%   cornea. For this reason, only first row of values are used by the
%   calling function.
%
% Inputs:
%   eyePoint
%   angle_p1p2, angle_p1p3 - Scalar in radians. The angle w.r.t. the optical
%                           axis of the initial ray. 
%   opticalSytem
%
% Outputs:
%   virtualEyePoint       - A 1x3 matrix that is the position of the
%                           virtual image.
%


% Ray trace for these thetas
outputRayEyeWorld = rayTraceEllipsoids(eyePoint, [angle_p1p2, angle_p1p3], opticalSystem);

% If we received a ray-trace error, then return nans for output ray
if isempty(outputRayEyeWorld) || isempty(outputRayEyeWorld)
    virtualEyePoint = nan(2,3);
    return
end

% Adjust the p1 (optical axis) position of the ray to have an initial
% position at the depth of the pupil
slope_p1p2 =(outputRayEyeWorld(2,2)-outputRayEyeWorld(1,2))/(outputRayEyeWorld(2,1)-outputRayEyeWorld(1,1));
slope_p1p3 =(outputRayEyeWorld(2,3)-outputRayEyeWorld(1,3))/(outputRayEyeWorld(2,1)-outputRayEyeWorld(1,1));
zOffset=outputRayEyeWorld(1,1)-eyePoint(1);
outputRayEyeWorld(:,1)=outputRayEyeWorld(:,1)-zOffset;
outputRayEyeWorld(:,2)=outputRayEyeWorld(:,2)-(zOffset*slope_p1p2);
outputRayEyeWorld(:,3)=outputRayEyeWorld(:,3)-(zOffset*slope_p1p3);

virtualEyePoint = outputRayEyeWorld(1,:);

end % calcVirtualImageRay


