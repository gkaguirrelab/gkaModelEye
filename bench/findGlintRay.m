function [outputRay, initialRay, targetIntersectError ] = findGlintRay( worldOrigin, eyePose, worldTarget, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR )
% Finds the ray that starts at worldOrigin and intersects worldTarget
%
% Syntax:
%  [outputRay, initialRay, targetIntersectError ] = findGlintRay( worldOrigin, eyePose, worldTarget, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR )
%
% Description:
%   This routine returns the outputRay from the last surface of an optical
%   system for a point that has originated from a worldOrigin and has
%   arrived at an observer positioned at the worldTarget location. The
%   routine accounts for rotation of the eye specified in the eyePose
%   variable. The outputRay returned by this routine provides the location
%   of the virtual image of that point.
%
%   The routine is designed for the particular case in which the
%   worldOrigin is an active light source associated with a camera located
%   at worldTarget, and the outputRay is reflected off the corneal surface,
%   creating a "glint" in the camera image.
%
%   the routine takes as input three optical systems. This allows a ray to
%   be traced from the light source 1) in the camera-to-eye direction
%   through any non-rotating refractive elements in the medium (e.g., a
%   spectacle lens), 2) reflecting off the rotating optical system (e.g.,
%   the tear film of the eye), 3) passing back through the non-rotating
%   elements, now in the eye-to-camera direction.
%
% Inputs:
%   worldOrigin           - A 3x1 vector that specifies the point in world 
%                           coordinates (x, y, z) of the source of a ray.
%                           This could be the coordinates of the light
%                           source of an active IR camera.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The stop
%                           radius value is unused by this routine.
%   worldTarget           - A 3x1 vector that specifies the point in world 
%                           coordinates (x, y, z) that the ray should
%                           intersect after exiting the optical system. A
%                           common application is to set worldTarget equal
%                           to the location of the pinhole aperture of a
%                           camera, which is found in:
%                           	sceneGeometry.cameraPosition.translation
%   rotationCenters       - Equal to sceneGeometry.eye.rotationCenters
%   opticalSystemFixRL    - Struct. This is the component of the optical
%                           system that is invariant with eye movements,
%                           traced in the camera-to-eye direction.
%                           Typically set to: sceneGeometry.
%                               refraction.cameraToMedium.opticalSystem
%   opticalSystemRot      - Struct. This is the component of the optical
%                           system that is subject to rotation with eye
%                           movements. Typically set to: sceneGeometry.
%                               refraction.stopToMedium.opticalSystem
%   opticalSystemFixLR    - Struct. This is the component of the optical
%                           system that is invariant with eye movements,
%                           traced in the eye-to-camera direction.
%                           Typically set to: sceneGeometry.
%                               refraction.mediumToCamera.opticalSystem
%
% Outputs:
%   outputRay             - 2x3 matrix that specifies in the eye coordinate
%                           space the ray as a unit vector of the form [p;
%                           u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity.
%                           Dimensions are p1, p2, p3.
%   initialRay            - A 2x3 vector that specifies in eye coordinate
%                           space the vector arising from the worldPoint
%                           that will intersect the worldTarget. Dimensions
%                           are p1, p2, p3.
%   targetIntersectError  - The distance (in mm) between the worldTarget
%                           and the closest passage of the outputRay
%
% Examples:
%{
    % Glint location in a rotated eye
    sceneGeometry = createSceneGeometry();
    eyePose = [-5, 3, 0, 2];
    cameraNodalPoint = sceneGeometry.cameraPosition.translation;
    irSourceLocation = cameraNodalPoint - sceneGeometry.cameraPosition.glintSourceRelative;
    rotationCenters = sceneGeometry.eye.rotationCenters;
    opticalSystemFixRL = sceneGeometry.refraction.cameraToMedium.opticalSystem;
    opticalSystemRot = sceneGeometry.refraction.glint.opticalSystem;
    opticalSystemFixLR = sceneGeometry.refraction.mediumToCamera.opticalSystem;
    [outputRay, initialRay, targetIntersectError ] = findGlintRay( irSourceLocation, eyePose, cameraNodalPoint, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR );
    outputRayCached = [ -0.090029371503151   1.140783928943272  -0.409912651321597;   0.994635325178541   0.088282493175748  -0.053914481409478 ];
    assert(max(max(abs(outputRay - outputRayCached))) < 1e-6)
%}


%% Variable and search option set-up

% Pre-define the output variables to keep the compiler happy
outputRay = nan(2,3);
initialRay = nan(2,3);
targetIntersectError = Inf;

% Set some parameters for the search
intersectErrorTolerance = 1e-4;
searchIterTolerance = 50;
searchIter = 0;
searchingFlag = true;
bestError = Inf;
bestAngles = [nan nan];

% Set fminbnd options to tolerate an error of 1e-2, and to make changes
% in theta as small as 1e-6.
TolFun = 1e-2; % intersection error to tolerate
TolX = 1e-6; % precision with which theta is estimated
options = optimset('TolFun',TolFun,'TolX',TolX,'Display','off');

% This structure will hold the rotation matrices to apply eye rotation. We
% define it here and then save the filled version of the variable that is
% returned by rotateEyeCoord. We won't need to compute it again, hopefully
% saving on execution time.
Rstruc = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);


%% Calculate initial guess for p1p2 and p1p3 angles
% Our goal is to find the angles with which a ray that departs worldOrigin
% results in the ray that ultimately intersects worldTarget. In the typical
% application, the worldOrigin is an active light source attached to a
% camera located at worldTarget, and the ray will have reflected off the
% cornea of the eye. Here we obtain an initial guess at the solution by
% finding the angles of a ray that connects the worldOrigin to the corneal
% apex

% When passed in the glint optical system, the cornea front surface is
% assigned a negative refractive index, which indicates that the ray-trace
% routine should reflect (instead of refract) the ray upon encountering
% this surface. We use this property to identify the reflective surface.
reflectSurfaceIdx = find(opticalSystemRot(:,end) < 0,1);

% Obtain the quadric matrix for this surface
S = quadric.vecToMatrix(opticalSystemRot(reflectSurfaceIdx,1:10));

% The apex of the ellipsoid has the geodetic coordinates [0 0]. We use this
% property to find the cartesian coordinates of the corneal apex. This
% solution holds even when the quadric surface for the cornea is slightly
% rotated out of alignment with the optical axis of the eye (i.e., if
% corneal tilt / tip rotations are being modeled). The location is
% specified in the eye coordinate space. The corneal apex can be located at
% either geodetic [0 0 0] or geodetic [0 180 0], depending upon the
% dimensions of the quadric. We test both and take the most anterior
% output.
cornealApexA = quadric.ellipsoidalGeoToCart([0 0 0],S)';
cornealApexB = quadric.ellipsoidalGeoToCart([0 180 0],S)';
if abs(cornealApexB(1)) < abs(cornealApexA(1))
    cornealApex = cornealApexB;
else
    cornealApex = cornealApexA;
end

% Next, we find the location of the corneal apex when the eye has been
% rotated by the angles specified in eyePose. We are still in the eye
% coordinate space.
[cornealApexRotated, Rstruc] = ...
    rotateEyeCoord(cornealApex, eyePose, rotationCenters, 'forward', Rstruc);

% Convert the origin point of the ray from world coordinates to eye
% coordinates.
eyeCoordOrigin = convertWorldToEyeCoord(worldOrigin);

% Find the angles with which a ray departs the ray origin and arrives at
% the corneal apex.
[angle_p1p2, angle_p1p3] = quadric.rayToAngles(quadric.normalizeRay([eyeCoordOrigin; cornealApexRotated - eyeCoordOrigin]'));

% If the absolute value of an initial search angle is greater than 90, we
% flip the direction in the search so that we don't get stuck at the
% -180/180 wrap around point.
p1p2Adjust = round(angle_p1p2/180)*180;
p1p3Adjust = round(angle_p1p3/180)*180;
angle_p1p2 = wrapTo180(angle_p1p2-p1p2Adjust);
angle_p1p3 = wrapTo180(angle_p1p3-p1p3Adjust);


%% Define the bounds for the search
% We use the angular extent of the bounding box for the reflective surface
% to set the upper and lower bounds on the possible solution angles.

% The bounding box for the reflective surface of the eye
bb = opticalSystemRot(reflectSurfaceIdx,12:17);

% Rotate the bounding box with the eye pose
bb([1 3 5])=rotateEyeCoord([bb(1) bb(3) bb(5)],eyePose, rotationCenters);
bb([2 4 6])=rotateEyeCoord([bb(2) bb(4) bb(6)],eyePose, rotationCenters);

% Find the bound angles by looping over the corners of the bounding box.
[~,p1Idx]=min(eyeCoordOrigin(1)-bb(1:2));
boundAngles_p1p2 = nan(2,2);
boundAngles_p1p3 = nan(2,2);
for xx=1:2
    for yy=1:2
        cornerPoint = [bb(p1Idx),bb(2+xx),bb(4+yy)];
        [p1p2B,p1p3B] = ...
            quadric.rayToAngles(quadric.normalizeRay([eyeCoordOrigin; cornerPoint-eyeCoordOrigin]'));
        boundAngles_p1p2(xx,yy) = wrapTo180(p1p2B-p1p2Adjust);
        boundAngles_p1p3(xx,yy) = wrapTo180(p1p3B-p1p3Adjust);
    end
end

% Take the max and mins of the angles to set bounds on the search
ub_p1p2 = max(max(boundAngles_p1p2));
ub_p1p3 = max(max(boundAngles_p1p3));
lb_p1p2 = min(min(boundAngles_p1p2));
lb_p1p3 = min(min(boundAngles_p1p3));


%% Create an anonymous function for ray tracing
intersectErrorFunc = @(p1p2,p1p3) calcTargetIntersectError(eyeCoordOrigin, wrapTo180(p1p2+p1p2Adjust), wrapTo180(p1p3+p1p3Adjust), eyePose, worldTarget, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc);


%% Shrink bounds
% Restrict to the domain of valid ray trace solutions
errorFunc = @(x) intersectErrorFunc(x,angle_p1p3);
ub_p1p2 = shrinkBound(angle_p1p2,ub_p1p2,TolX,errorFunc);
lb_p1p2 = shrinkBound(angle_p1p2,lb_p1p2,TolX,errorFunc);
errorFunc = @(x) intersectErrorFunc(angle_p1p2,x);
ub_p1p3 = shrinkBound(angle_p1p3,ub_p1p3,TolX,errorFunc);
lb_p1p3 = shrinkBound(angle_p1p3,lb_p1p3,TolX,errorFunc);


%% Perform the search

% Get the intial target error
targetIntersectError = intersectErrorFunc(angle_p1p2, angle_p1p3);

% If the x0 guess isn't good enough, proceed with the search
if targetIntersectError > intersectErrorTolerance
    
    % Enter a while loop that iteratively refines the theta values until
    % criteria are met.
    while searchingFlag

        % Update the last error
        lastError = targetIntersectError;
        
        % The distance error function for searching across p1p2 values
        myObj_p1p2 = @(angleX) intersectErrorFunc(angleX, angle_p1p3);
        
        % Perform the search across angle_p1p2
        [angle_p1p2,~]=fminbnd(myObj_p1p2,lb_p1p2,ub_p1p2,options);
                
        % We hit a bound on the first iteration. Switch to a grid search
        % for x0 and then fminsearch.
        if searchIter == 0 && any([...
                abs(angle_p1p2 - lb_p1p2)<1e-3,...
                abs(angle_p1p2 - ub_p1p2)<1e-3])
            [xxVals,yyVals]=meshgrid(-75:5:75,-75:5:75);
            gridError = inf(length(xxVals),length(xxVals));
            for xx=1:length(xxVals)
                for yy=1:length(yyVals)
                    gridError(xx,yy)=intersectErrorFunc(xxVals(xx,yy),yyVals(xx,yy));
                end
            end
            [~,idx]=min(gridError(:));
            angle_p1p2 = xxVals(idx);
            angle_p1p3 = yyVals(idx);
            myObj = @(p) intersectErrorFunc(p(1),p(2));
            [p,targetIntersectError]=fminsearch(myObj,[angle_p1p2 angle_p1p3]);
            angle_p1p2 = p(1); angle_p1p3 = p(2);
            searchingFlag = false;
            continue
        end
        
        % The distance error function for searching across p1p3 values
        myObj_p1p3 = @(angleX) intersectErrorFunc(angle_p1p2, angleX);
        
        % Perform the search across angle_p1p3
        [angle_p1p3,targetIntersectError]=fminbnd(myObj_p1p3,lb_p1p3,ub_p1p3,options);

        % Update the best result
        if targetIntersectError < bestError
            bestError = targetIntersectError;
            bestAngles = [angle_p1p2 angle_p1p3];
        end
        
        % Iterate the search count
        searchIter = searchIter+1;
                
        % Determine if we have met a stopping criterion
        if targetIntersectError<intersectErrorTolerance || (lastError-targetIntersectError)<intersectErrorTolerance || searchIter>searchIterTolerance
            angle_p1p2 = bestAngles(1); angle_p1p3 = bestAngles(2);
            targetIntersectError = bestError;
            searchingFlag = false;            
        end
        
    end % while searching
    
end % Test x0 guess


% Remove any angle adjustment
angle_p1p2 = wrapTo180(angle_p1p2+p1p2Adjust);
angle_p1p3 = wrapTo180(angle_p1p3+p1p3Adjust);


%% Obtain the initial and output rays

% Assemble the input ray. Note that the rayTraceQuadrics routine handles
% vectors as a 3x2 matrix, as opposed to a 2x3 matrix in this function.
% Tranpose operations ahead.
initialRay = quadric.anglesToRay(eyeCoordOrigin', angle_p1p2, angle_p1p3)';

% Trace the system using the initialRay defined by the solution angles
outputRay = threeSystemTrace(initialRay', eyePose, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc);

% Counter-rotate the output ray to account for eye rotation
outputRay = rotateEyeRay(outputRay, eyePose, rotationCenters, 'inverse');


end % findGlintRay -- MAIN




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Angle wrap functions
% Shadowing the Mathworks code to allow code generation
function lon = wrapTo180(lon)
q = (lon < -180) | (180 < lon);
lon(q) = wrapTo360(lon(q) + 180) - 180;
end

function lon = wrapTo360(lon)
positiveInput = (lon > 0);
lon = mod(lon, 360);
lon((lon == 0) & positiveInput) = 360;
end


%% shrinkBound
% Use divide approach to find a bound that returns a value for the passed
% function that is less than 1e6
function boundOut = shrinkBound(x0,boundIn,TolX,errorFunc)

% Check if we already have a valid bound, in which case return
if errorFunc(boundIn) < 1e6
    boundOut = boundIn;
    return
end

% Set up our vars
step = -(boundIn-x0)/2;
boundOut = boundIn + step;
stillSearching = true;
lastValidVal = x0;
lastBound = boundIn;

% Let's go searching
while stillSearching
    if errorFunc(boundOut) < 1e6
        % We have a valid trace.
        if abs(step) < TolX
            % If our step is below the TolX, then we are
            % done.
            stillSearching = false;
        else
            % We are valid, but still too coarse, so move half-way back
            % towards the lastBound
            lastValidVal = boundOut;
            step = +(lastBound-boundOut)/2;
            lastBound = boundOut;
            boundOut = boundOut + step;
        end
    else
        % Not a valid trace
        if abs(step) < TolX
            % We are below the search tol, so report the last valid bound
            % we found
            boundOut = lastValidVal;
            stillSearching = false;
        else
            % Move half-way towards the x0
            step = -(boundOut-x0)/2;
            lastBound = boundOut;
            boundOut = boundOut + step;
        end
    end
end
end


%% calcTargetIntersectError
function distance = calcTargetIntersectError(eyeCoordOrigin, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc)
% Smallest distance of the exit ray from worldTarget
%
% Syntax:
%  distance = calcTargetIntersectError(eyePoint, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystem)
%
% Description:
%   This function returns the Euclidean distance between a target in the
%   world coordinate system and a ray that has exited from the optical
%   system of a rotated eye. This distance is calculated within eyeWorld
%   coordinates.
%
%   This function is used to find angles in the p1p2 and p1p3 planes that
%   minimize the distance between the intersection point of the ray in the
%   camera plane and the pinhole aperture of a camera. At a distance of
%   zero, the ray would enter the pinhole aperture of the camera and thus
%   produce a point on the resulting image.
%
% Inputs:
%   eyeCoordOrigin        - The origin of the ray in eye coordinates
%   angle_p1p2, angle_p1p3 - Scalars in degrees. The angle w.r.t. the 
%                           optical axis of the initial ray. 
%   eyePose               - As defined in the main function.
%   worldTarget
%   rotationCenters
%   opticalSystem
%
% Outputs:
%   distance              - Scalar in units of mm. The minimum Euclidean 
%                           distance between the worldTarget and a ray
%                           exiting the optical system of the rotated eye.
%                           Set to Inf if an error is returned by
%                           rayTraceQuadrics.
%


% Assemble the input ray. Note that the rayTraceQuadrics routine handles
% vectors as a 3x2 matrix, as opposed to a 2x3 matrix in this function.
% Tranpose operations ahead.
inputRayEyeCoord = quadric.anglesToRay(eyeCoordOrigin',angle_p1p2,angle_p1p3);

% Conduct the ray trace through the fixed, rotating, and back through the
% fixed optical systems
outputRayEyeCoord = threeSystemTrace(inputRayEyeCoord, eyePose, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc);

% If any must intersect surfaces were missed, the output ray will contain
% nans. In this case, return 1e6 for the distance
if any(isnan(outputRayEyeCoord))
    distance = 1e6;
    return
end

% Put the target point in eye coordinate space.
eyeCoordTarget = convertWorldToEyeCoord(worldTarget);

% Calculate the distance between the closest approach of the outputRay to
% the target.
distance = quadric.distancePointRay(eyeCoordTarget',outputRayEyeCoord');

end % calcTargetIntersectError



function outputRayEyeCoord = threeSystemTrace(inputRayEyeCoord, eyePose, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc)

% Ray trace through the fixed system (typically, camera to the medium
% adjacent to the eye)
outputRayEyeCoord = rayTraceQuadrics(inputRayEyeCoord, opticalSystemFixRL)';

% Counter-rotate the outputRayEyeCoord by the eye pose. This places the
% outputRayEyeCoord so that the eyeCoordTarget is in a position equivalent
% to if the eye had rotated.
outputRayEyeCoord = rotateEyeRay(outputRayEyeCoord, eyePose, rotationCenters, 'inverse', Rstruc);

% Ray trace through the eye system that is subject to rotation
outputRayEyeCoord = rayTraceQuadrics(outputRayEyeCoord', opticalSystemRot)';

% Return the ray to the original eyeWorld coordinate frame, which is
% aligned with the optical axis of the eye
outputRayEyeCoord = rotateEyeRay(outputRayEyeCoord, eyePose, rotationCenters, 'forward', Rstruc);

% Ray trace back through the fixed system that is not subject to rotation
outputRayEyeCoord = rayTraceQuadrics(outputRayEyeCoord', opticalSystemFixLR)';


end
