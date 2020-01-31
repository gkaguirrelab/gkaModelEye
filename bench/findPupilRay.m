function [outputRay, initialRay, targetIntersectError ] = findPupilRay( eyePoint, eyePose, worldTarget, rotationCenters, opticalSystem )
% Returns the virtual image ray of a point in eyeWorld coordinates
%
% Syntax:
%  [outputRay, initialRay, targetIntersectError ] = findPupilRay( eyePoint, eyePose, worldTarget, rotationCenters, opticalSystem )
%
% Description:
%   This routine returns the outputRay from the last surface of an optical
%   system for a point that has originated from an eyeWorld coordinate
%   point and has arrived at an observer positioned at the worldTarget
%   location. The routine accounts for rotation of the eye specified in the
%   eyePose variable. The outputRay returned by this routine provides the
%   location of the virtual image of that point.
%
% Inputs:
%   eyePoint              - A 1x3 vector that gives the coordinates (in mm)
%                           of a point in eyeWorld space with the
%                           dimensions p1, p2, p3.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees. The stop
%                           radius value is unused by this routine.
%   worldTarget     -       A 3x1 vector that specifies the point in world 
%                           coordinates (x, y, z) that the ray should
%                           intersect after exiting the optical system. A
%                           common application is to set worldTarget equal
%                           to the location of the pinhole aperture of a
%                           camera, which is found in:
%                           	sceneGeometry.cameraPosition.translation
%   rotationCenters       - Equal to sceneGeometry.eye.rotationCenters
%   opticalSystem         - Typically set equal to: sceneGeometry.
%                               refraction.stopToCamera.opticalSystem
%
% Outputs:
%   outputRay             - 2x3 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%                           dimensions p1, p2, p3.
%   initialRay            - A 2x3 vector that specifies in eyeWorld space
%                           the vector arising from the eyePoint that will
%                           intersect the worldTarget.
%   targetIntersectError  - The distance (in mm) between the worldTarget
%                           and the closest passage of a ray arising from
%                           the eyeWorld point after it exits the optical
%                           system.
%
% Examples:
%{
    % Basic example that finds the virtual image location for the center of
    % the aperture stop, with eye and the camera in their default positions
    sceneGeometry = createSceneGeometry();
    % Assemble the args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
    	sceneGeometry.eye.rotationCenters, ...
    	sceneGeometry.refraction.stopToCamera.opticalSystem};
    outputRay = findPupilRay( sceneGeometry.eye.stop.center, [-5 10 0 2], args{:} );
    % Test output against cached value
    outputRayCached = [ -0.028495571780279   0.316608707224537  -0.630653051462113; 0.978179527188788   0.093602909096939  -0.185481287458627];
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
    	[~,~,~,~,~,pointLabels,errors]=pupilProjection_fwd(eyePoses(pp,:),sceneGeometry,'nStopPerimPoints',nStopPerimPoints,'rayTraceErrorThreshold',rayTraceErrorThreshold);
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



%% Find the p1p2 and p1p3 angles
% For this eyeWorld point, we find the angles of origin of a ray in the
% p1p2 and p1p3 planes that results in an exit ray (after passing through
% the optical system) that passes as close as possible to the worldTarget

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

% Set the inital guess for the angles by finding (w.r.t. the optical axis)
% the angle of the ray that connects the eye point to the worldTarget
% (after re-arranging the dimensions of the worldTarget variable).
eyeCoordTarget = relativeCameraPosition(eyePose, worldTarget, rotationCenters);
[angle_p1p2, angle_p1p3] = quadric.rayToAngles(quadric.normalizeRay([eyePoint; eyeCoordTarget-eyePoint]'));

% Set bounds on the search
ub_p1p2 = 90;
ub_p1p3 = 90;
lb_p1p2 = -90;
lb_p1p3 = -90;

% Ensure that the initial guess is in bounds
angle_p1p2 = max([angle_p1p2 lb_p1p2]);
angle_p1p2 = min([angle_p1p2 ub_p1p2]);
angle_p1p3 = max([angle_p1p3 lb_p1p3]);
angle_p1p3 = min([angle_p1p3 ub_p1p3]);

% Create an anonymous function for ray tracing
intersectErrorFunc = @(p1p2,p1p3) calcTargetIntersectError(eyePoint, p1p2, p1p3, eyePose, worldTarget, rotationCenters, opticalSystem);

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
initialRay = quadric.anglesToRay(eyePoint', angle_p1p2, angle_p1p3)';

% Ray trace
outputRay = rayTraceQuadrics(initialRay', opticalSystem);
outputRay = outputRay';

end % findPupilRay -- MAIN




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calcTargetIntersectError
function distance = calcTargetIntersectError(eyePoint, angle_p1p2, angle_p1p3, eyePose, worldTarget, rotationCenters, opticalSystem)
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
%   eyePoint
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
inputRay = quadric.anglesToRay(eyePoint',angle_p1p2,angle_p1p3);

% Ray trace
outputRayEyeWorld = rayTraceQuadrics(inputRay, opticalSystem);
outputRayEyeWorld = outputRayEyeWorld';

% If any must intersect surfaces were missed, the output ray will contain
% nans. In this case, return Inf for the distance
if any(isnan(outputRayEyeWorld))
    distance = Inf;
    return
end

% We assign the variable eyeCoordTarget the coordinates of the target after
% conversion to eye coordinates. Then, the point is counter-rotated by the
% eye pose, so that the eyeCoordTarget is in a position equivalent to if
% the eye had rotated.
eyeCoordTarget = relativeCameraPosition(eyePose, worldTarget, rotationCenters);

% Calculate the distance between the closest approach of the outputRay to
% the target.
distance = quadric.distancePointRay(eyeCoordTarget',outputRayEyeWorld');

end % calcTargetIntersectError



function eyeCoordTarget = relativeCameraPosition(eyePose, worldTarget, rotationCenters)

% We assign the variable eyeCoordTarget the coordinates of the target after
% conversion to eye coordinates. Then, the point is counter-rotated by the
% eye pose, so that the eyeCoordTarget is in a position equivalent to if
% the eye had rotated.
cameraRot = -eyePose;
RotAzi = [cosd(cameraRot(1)) -sind(cameraRot(1)) 0; sind(cameraRot(1)) cosd(cameraRot(1)) 0; 0 0 1];
RotEle = [cosd(-cameraRot(2)) 0 sind(-cameraRot(2)); 0 1 0; -sind(-cameraRot(2)) 0 cosd(-cameraRot(2))];
RotTor = [1 0 0; 0 cosd(cameraRot(3)) -sind(cameraRot(3)); 0 sind(cameraRot(3)) cosd(cameraRot(3))];

% Rearrange the worldTarget dimensions to switch from world to eye
% coordinate space. This is now the eyeTarget or ET
eyeCoordTarget = worldTarget([3 1 2])';

% Torsion
eyeCoordTarget=eyeCoordTarget-rotationCenters.tor;
eyeCoordTarget = (RotTor*eyeCoordTarget')';
eyeCoordTarget=eyeCoordTarget+rotationCenters.tor;
% Elevation
eyeCoordTarget=eyeCoordTarget-rotationCenters.ele;
eyeCoordTarget = (RotEle*eyeCoordTarget')';
eyeCoordTarget=eyeCoordTarget+rotationCenters.ele;
% Azimuth
eyeCoordTarget=eyeCoordTarget-rotationCenters.azi;
eyeCoordTarget = (RotAzi*eyeCoordTarget')';
eyeCoordTarget=eyeCoordTarget+rotationCenters.azi;
end
