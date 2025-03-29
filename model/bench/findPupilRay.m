function [outputRay, initialRay, targetIntersectError ] = findPupilRay( eyeCoordOrigin, eyePose, worldTarget, rotationCenters, opticalSystemRot, opticalSystemFix)
% Finds the ray that starts at eyeCoordOrigin and intersects worldTarget
%
% Syntax:
%  [outputRay, initialRay, targetIntersectError ] = findPupilRay( eyeCoordOrigin, eyePose, worldTarget, rotationCenters, opticalSystemRot, opticalSystemFix )
%
% Description:
%   This routine returns the outputRay from the last surface of an optical
%   system for a point that has originated from an eyeWorld coordinate and
%   has arrived at an observer positioned at the worldTarget location. The
%   routine accounts for rotation of the eye specified in the eyePose
%   variable. The outputRay returned by this routine provides the location
%   of the virtual image of that point.
%
% Inputs:
%   eyeCoordOrigin        - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in the eye coordinate space with the
%                           dimensions p1, p2, p3.
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
%   opticalSystemRot      - Struct. This is the component of the optical
%                           system that is subject to rotation with eye
%                           movements. Typically set to: sceneGeometry.
%                               refraction.stopToMedium.opticalSystem
%   opticalSystemFix      - Struct. This is the component of the optical
%                           system that is invariant with eye movements.
%                           Typically set to: sceneGeometry.
%                               refraction.mediumToCamera.opticalSystem
%
% Outputs:
%   outputRay             - 2x3 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%                           Dimensions p1, p2, p3.
%   initialRay            - A 2x3 vector that specifies in eyeWorld space
%                           the vector arising from the eyeCoordOrigin that will
%                           intersect the worldTarget.
%   targetIntersectError  - The distance (in mm) between the worldTarget
%                           and the closest passage of the outputRay.
%
% Examples:
%{
    % Basic example that finds the virtual image location for the center of
    % the aperture stop, with eye and the camera in their default positions
    sceneGeometry = createSceneGeometry();
    % Assemble the args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.stopToMedium.opticalSystem, ...
        sceneGeometry.refraction.mediumToCamera.opticalSystem};
    outputRay = findPupilRay( sceneGeometry.eye.stop.center, [-5 10 0 2], args{:} );
    % Test output against cached value
    outputRayCached = [ -0.027365339254198   0.320630761081699  -0.625739056740334;  0.978152533979250   0.095064502933833  -0.184879313477513 ];
    assert(max(max(abs(outputRay - outputRayCached))) < 1e-6)
%}
%{
    %% Confirm that targetIntersectError remains small across eye poses
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Perform 100 forward projections with randomly selected eye poses
    nPoses = 100;
    eyePoses=[(rand(nPoses,1)-0.5)*70, (rand(nPoses,1)-0.5)*70, zeros(nPoses,1), 2+(rand(nPoses,1)-0.5)*1];
    nStopPerimPoints = 6;
    rayTraceErrorThreshold = 1;
    clear targetIntersectError
    targetIntersectError = nan(nStopPerimPoints,nPoses);
    for pp = 1:nPoses
    	[~,~,~,~,~,~,pointLabels,errors]=projectModelEye(eyePoses(pp,:),sceneGeometry,'nStopPerimPoints',nStopPerimPoints,'rayTraceErrorThreshold',rayTraceErrorThreshold);
        idx = strcmp(pointLabels,'pupilPerimeter');
        targetIntersectError(1:sum(idx),pp) = errors(idx);
    end
    % Make sure the targetIntersectError is small and not systematically
    % related to eyePose
    figure
    plot(sqrt(eyePoses(:,1).^2+eyePoses(:,2).^2),nanmax(targetIntersectError),'.r')
    xlabel('Euclidean rotation distance [deg]');
    ylabel('Max ray trace camera intersection error [mm]');
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
% returned by rotateEyeCoord. We then don't need to compute it again,
% hopefully saving on execution time.
Rstruc = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);

% Pre-covert the optical system surfaces into quadric matrix form
tensorSRot = nan(4,4,100);
nSurfacesRot = sum(~all(isnan(opticalSystemRot),2));
for ii = 1:nSurfacesRot
    tensorSRot(:,:,ii) = quadric.vecToMatrix(opticalSystemRot(ii,1:10));
end
tensorSFix = nan(4,4,100);
nSurfacesFix = sum(~all(isnan(opticalSystemFix),2));
for ii = 1:nSurfacesFix
    tensorSFix(:,:,ii) = quadric.vecToMatrix(opticalSystemFix(ii,1:10));
end

%% Calculate initial guess for p1p2 and p1p3 angles
% Our goal is to find the angles with which a ray departs eyeCoordOrigin
% that results in a ray that ultimately intersects worldTarget. In the
% typical application, eyeCoordOrigin is a point on the border of the
% aperture stop of the iris, and  worldTarget is the nodal point of a
% camera that is observing the eye.  Here we obtain an initial guess at the
% solution by finding the angles of a ray that connects the eyeCoordOrigin
% to the worldTarget


% Convert the target from world to eye coordinates
eyeCoordTarget = convertWorldToEyeCoord(worldTarget);

% We subject the target to an inverse rotation with the angles specified by
% the eyePose. That is, the location of the camera relative to the eye is
% adjusted to account for the rotation of the eye
[eyeCoordTargetRotated, Rstruc] = rotateEyeCoord(eyeCoordTarget, eyePose, rotationCenters, 'inverse', Rstruc);

% Find the angles that connect a ray from the origin to the target
[angle_p1p2, angle_p1p3] = quadric.rayToAngles(quadric.normalizeRay([eyeCoordOrigin; eyeCoordTargetRotated-eyeCoordOrigin]'));


%% Define the bounds for the search

% Set bounds on the search
ub_p1p2 = 90;
ub_p1p3 = 90;
lb_p1p2 = -90;
lb_p1p3 = -90;

% Create an anonymous function for ray tracing
intersectErrorFunc = @(p1p2,p1p3) calcTargetIntersectError(...
    eyeCoordOrigin, p1p2, p1p3, eyePose, worldTarget, rotationCenters, ...
    opticalSystemRot, opticalSystemFix, ...
    tensorSRot, tensorSFix, ...
    nSurfacesRot, nSurfacesFix, ...
    Rstruc);

% Shrink the bounds to restrict to the domain of valid ray trace solutions
errorFunc = @(x) intersectErrorFunc(x,angle_p1p3);
ub_p1p2 = shrinkBound(angle_p1p2,ub_p1p2,TolX,errorFunc);
lb_p1p2 = shrinkBound(angle_p1p2,lb_p1p2,TolX,errorFunc);
errorFunc = @(x) intersectErrorFunc(angle_p1p2,x);
ub_p1p3 = shrinkBound(angle_p1p3,ub_p1p3,TolX,errorFunc);
lb_p1p3 = shrinkBound(angle_p1p3,lb_p1p3,TolX,errorFunc);

% Ensure that the initial guess is in bounds
angle_p1p2 = max([angle_p1p2 lb_p1p2]);
angle_p1p2 = min([angle_p1p2 ub_p1p2]);
angle_p1p3 = max([angle_p1p3 lb_p1p3]);
angle_p1p3 = min([angle_p1p3 ub_p1p3]);


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


%% Obtain the initial and output rays

% Assemble the input ray. Note that the rayTraceQuadrics routine handles
% vectors as a 3x2 matrix, as opposed to a 2x3 matrix in this function.
% Tranpose operations ahead.
initialRay = quadric.anglesToRay(eyeCoordOrigin', angle_p1p2, angle_p1p3);

% Trace the system using the initialRay defined by the solution angles
outputRay = twoSystemTrace(initialRay, eyePose, rotationCenters, ...
    opticalSystemRot, opticalSystemFix, ...
    tensorSRot, tensorSFix, ...
    nSurfacesRot, nSurfacesFix, ...
    Rstruc);

% Undo the effect of eye rotation to obtain the outputRay in the coordinate
% space of the eye prior to rotation
outputRay = rotateEyeRay(outputRay, eyePose, rotationCenters, 'inverse');

% Transpose the initialRay to return
initialRay = initialRay';

end % findPupilRay -- MAIN




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
function distance = calcTargetIntersectError(eyeCoordOrigin, ...
    angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, ...
    opticalSystemRot, opticalSystemFix, ...
    tensorSRot, tensorSFix, ...
    nSurfacesRot, nSurfacesFix, ...
    Rstruc)
% Smallest distance of the exit ray from worldTarget
%
% Syntax:
%  distance = calcTargetIntersectError(eyeCoordOrigin, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystemRot, opticalSystemFix, tensorSRot, tensorSFix, Rstruc)
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
%   eyeCoordOrigin
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

% Test the input angles. If they are greater than +-90, return inf
if abs(angle_p1p2)>=90 || abs(angle_p1p3)>=90
    distance = Inf;
    return
end

% Assemble the input ray. Note that the rayTraceQuadrics routine handles
% vectors as a 3x2 matrix, as opposed to a 2x3 matrix in this function.
% Tranpose operations ahead.
inputRayEyeCoord = quadric.anglesToRay(eyeCoordOrigin',angle_p1p2,angle_p1p3);

% Conduct the ray trace through the rotating and fixed optical systems
outputRayEyeCoord = twoSystemTrace(inputRayEyeCoord, eyePose, rotationCenters, ...
    opticalSystemRot, opticalSystemFix, ...
    tensorSRot, tensorSFix, ...
    nSurfacesRot, nSurfacesFix, ...
    Rstruc);

% If any must intersect surfaces were missed, the output ray will contain
% nans. In this case, return Inf for the distance
if any(isnan(outputRayEyeCoord))
    distance = Inf;
    return
end

% Move the worldTarget into the eye coordinate space
eyeCoordTarget = convertWorldToEyeCoord(worldTarget);

% Calculate the distance between the closest approach of the outputRay to
% the target.
distance = quadric.distancePointRay(eyeCoordTarget',outputRayEyeCoord');

end % calcTargetIntersectError


function outputRayEyeCoord = twoSystemTrace(inputRayEyeCoord, eyePose, rotationCenters, ...
    opticalSystemRot, opticalSystemFix, ...
    tensorSRot, tensorSFix, ...
    nSurfacesRot, nSurfacesFix, ...
    Rstruc)

% Ray trace through the eye system that is subject to rotation
outputRayEyeCoord = rayTraceQuadrics(inputRayEyeCoord, opticalSystemRot, tensorSRot, nSurfacesRot)';

% Apply the eye rotation.
outputRayEyeCoord = rotateEyeRay(outputRayEyeCoord, eyePose, rotationCenters, 'forward', Rstruc);

% Ray trace through the fixed system that is not subject to rotation
outputRayEyeCoord = rayTraceQuadrics(outputRayEyeCoord', opticalSystemFix, tensorSFix, nSurfacesFix)';


end


