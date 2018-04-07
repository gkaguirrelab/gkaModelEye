function [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, opticalSystem, extrinsicTranslationVector, rotationCenters )
% Returns the virtual image location of a point in eyeWorld coordinates
%
% Syntax:
%  [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( eyeWorldPoint, eyePose, opticalSystem, extrinsicTranslationVector, rotationCenters )
%
% Description:
%   This routine returns a handle to a function that is used to calculate
%   the location of a virtual image point that has passed through an
%   optical system. If the key-value 'functionDirPath' is set, then the
%   handle will be to a compiled mex function that has been written to
%   disk, otherwise, the handle will be to a matlab function in memory.
%   The former executes ~100x faster.
%
%   The virtualImageFunc is assembled step-wise from more elementary
%   algorithms:
%       traceOpticalSystem - 2D ray tracing through the cornea and any
%           corrective lenses
%       calcCameraNodeDistanceError2D - 2D distance of ray intersection on 
%           camera plane from camera node
%       calcVirtualImageRay - Returns the unit vector virtual image ray for
%           the initial depth position
%   This main routine calls out to local functions that assemble these
%   elementary components. Each elementary component is either saved as a
%   file at the specified location, or maintained as a function in memory.
%   Once these components are assembled, a handle is made to the function
%   virtualImageFunc, which uses these elementary components.
%
% Inputs:
%   sceneGeometry         - A sceneGeometry structure. Critically, this
%                           includes an optical system.
%
% Optional key-value pairs:
%  'functionDirPath'     - Character string, default empty. If set this
%                           defines the location in which the compiled
%                           function is writen.
%  'cleanUpCompileDir'    - Logical, default true. If file path is
%                           provided, this flag determines if the
%                           intermediate compilation products are deleted.
%
% Outputs:
%   virtualImageFuncPointer - Structure. Includes the fields:
%                           'handle' - handle for the function.
%                           'path' -  full path to the stored mex file; set
%                               to empty if stored only in memory.
%                           'opticalSystem' - the optical system used to
%                               generate the function.
%
% Examples:
%{
    % Basic example, placing the function in memory
    sceneGeometry = createSceneGeometry();
    [virtualEyeWorldPoint, nodalPointIntersectError] = virtualImageFunc( [-3.7, 2 0], [0 0 0], sceneGeometry.opticalSystem, sceneGeometry.extrinsicTranslationVector, sceneGeometry.eye.rotationCenters )

%}
%{
    % Basic example with file caching of the functions
    sceneGeometry = createSceneGeometry();
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry, 'functionDirPath', '/tmp/demo_virtualImageFunc' );
%}
%{
    % Demonstrate how the time it takes to perform the symbolic variable
    % calculations grows geometrically with the number of surfaces in the
    % optical system.

    % Obtain a default sceneGeometry. 
    sceneGeometry = createSceneGeometry();
    % Define the virtual image function
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(1)=toc;
    n(1)=size(sceneGeometry.opticalSystem,1);
    % Add a contact lens (one additional surface)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
    % Define the ray tracing functions 
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(2)=toc;
    n(2)=size(sceneGeometry.opticalSystem,1);
    % Add a spectacle lens (two additional surfaces)
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',-2);
    % Define the ray tracing functions 
    tic
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
    t(3)=toc;
    n(3)=size(sceneGeometry.opticalSystem,1);
    % Plot the timing results
    plot(n,t,'*r');
    xlabel('# of surfaces in optical model');
    ylabel('time to assemble ray tracing funcs [secs]');
%}
% 
%     
%      traceOpticalSystem = matlabFunction(outputRayFunc, ...
%          'Vars',{z h theta});
    

traceOpticalSystem = @(z, h, theta) rayTraceCenteredSurfaces([z h], theta,opticalSystem);
cameraNodeDistanceError2D_p1p2 = @(theta_p1p2) calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, extrinsicTranslationVector, rotationCenters, traceOpticalSystem);
cameraNodeDistanceError3D = @(theta_p1p2, theta_p1p3) calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, extrinsicTranslationVector, rotationCenters, traceOpticalSystem);
virtualImageRay = @(theta_p1p2, theta_p1p3) calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, traceOpticalSystem);


options = optimset('TolFun',1e-2,'TolX',1e-6);

% Define an error function which is the distance between the nodal
% point of the camera and the point at which a ray intersects the
% plane that contains the camera, with the ray departing from the
% eyeWorld point at angle theta in the p1p2 plane.
errorFunc = @(theta_p1p2) cameraNodeDistanceError2D_p1p2(theta_p1p2);
% Conduct an fminsearch to find the p1p2 theta that results in a
% ray that strikes as close as possible to the camera nodal point.
% Because the errorFunc returns nan for values very close to zero,
% we initialize the search with a slightly non-zero value (1e-4)
theta_p1p2=fminsearch(errorFunc,1e-4,options);
% Now repeat this process for in the p1p3 plane, setting the p1p2 plane to
% the theta value that was just found
errorFunc = @(theta_p1p3) cameraNodeDistanceError3D(theta_p1p2, theta_p1p3);
% The fVal at the solution is the the total error (in mm) in both
% dimensions for intersecting the nodal point of the camera.
[theta_p1p3, nodalPointIntersectError]=fminsearch(errorFunc,1e-4,options);
% With both theta values calculated, now obtain the virtual image
% ray arising from the pupil plane that reflects the corneal optics
virtualImageRay = virtualImageRay(theta_p1p2, theta_p1p3);
% Extract the origin of the ray, which is the virtual image eyeWorld point
virtualEyeWorldPoint = virtualImageRay(1,:);



end % compileVirtualImageFunc -- MAIN



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% calcCameraNodeDistanceError2D_p1p2
function distance = calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, extrinsicTranslationVector, rotationCenters, traceOpticalSystem)
% 2D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = calcCameraNodeDistanceError2D_p1p2(eyeWorldPoint, theta_p1p2, eyePose, sceneGeometry, traceOpticalSystem)
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
%   We perform searches separately for thetas in the p1p2 and p1p3 planes
%   to minimize distance. This allows us to use fminsearch over a single
%   variable, which is computationally efficient.
%
% Inputs:
%   eyeWorldPoint
%   theta_p1p2            - Scalar in units of radians. The angle WRT the
%                           optical axis of the initial ray. The function
%                           is undefined for theta = 0 (i.e., a paraxial
%                           ray) and will return nan. Also, absolute values
%                           of pi correspond to a vertical ray that would
%                           not intersect with the optical system and thus
%                           will return nan. There are other combinations
%                           of eyeWorld positions and thetas that will
%                           return nan given the particular path of the ray
%                           through the optical system.
%   eyePose
%   rotationCenters
%   traceOpticalSystem
%
%
% Outputs:
%   distance              - Scalar in units of mm. The Euclidean distance
%                           of the intersection point of the ray on the
%                           Z camera plane from the nodal point of the
%                           camera.
%


outputRayEyeWorld2D_p1p2 = traceOpticalSystem(eyeWorldPoint(1), eyeWorldPoint(2),theta_p1p2);

% Add the third, constant dimension for the output rays
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

%% Torsion
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


%% Elevation
% p1p2
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


%% Azimuth
% p1p2
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

% Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
% point at which the ray will strike the plane that contains the camera
slope_xZ =(outputRaySceneWorld_p1p2(2,1)-outputRaySceneWorld_p1p2(1,1))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
slope_yZ =(outputRaySceneWorld_p1p2(2,2)-outputRaySceneWorld_p1p2(1,2))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
cameraPlaneX = outputRaySceneWorld_p1p2(1,1)+((extrinsicTranslationVector(3)-outputRaySceneWorld_p1p2(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p2(1,2)+((extrinsicTranslationVector(3)-outputRaySceneWorld_p1p2(1,3))*slope_yZ);

distance = sqrt((extrinsicTranslationVector(1)-cameraPlaneX)^2 + ...
        (extrinsicTranslationVector(2)-cameraPlaneY)^2 );


end % calcCameraNodeDistanceError2D



%% calcCameraNodeDistanceError3D
function distance = calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, extrinsicTranslationVector, rotationCenters, traceOpticalSystem)
% 3D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = calcCameraNodeDistanceError3D(eyeWorldPoint, theta_p1p2, theta_p1p3, eyePose, sceneGeometry, traceOpticalSystem)
%
% Description:
%   This function is similar to the calcCameraNodeDistanceError2D_p1p2,
%   except that it takes as input theta in both the p1p2 and p1p3 planes.
%   The distance value that is returned is still the Euclidean distance
%   between the intersection point of the output ray on the Z camera plane
%   and the nodal point of the camera.


outputRayEyeWorld2D_p1p2 = traceOpticalSystem(eyeWorldPoint(1), eyeWorldPoint(2), theta_p1p2);
outputRayEyeWorld2D_p1p3 = traceOpticalSystem(eyeWorldPoint(1), eyeWorldPoint(3), theta_p1p3);

% Create a 3D output ray system. Shift the p1p3 ray to have the same
% initial p1 value as the p1p2 vector
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
% Prepare to rotate the outputRay into the sceneWorld coordinates
RotAzi = [cosd(eyePose(1)) -sind(eyePose(1)) 0; sind(eyePose(1)) cosd(eyePose(1)) 0; 0 0 1];
RotEle = [cosd(eyePose(2)) 0 sind(eyePose(2)); 0 1 0; -sind(eyePose(2)) 0 cosd(eyePose(2))];
RotTor = [1 0 0; 0 cosd(eyePose(3)) -sind(eyePose(3)); 0 sind(eyePose(3)) cosd(eyePose(3))];

% Copy over the outputRay from eye to head world
outputRayHeadWorld3D=outputRayEyeWorld3D;

%% Torsion
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


%% Elevation
% p1p2
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


%% Azimuth
% p1p2
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

% We reverse the direction of the Y axis so that positive elevation of the
% eye corresponds to a movement of the pupil upward in the image
outputRaySceneWorld3D(:,2) = outputRaySceneWorld3D(:,2)*(-1);

% Obtain an expression for X and Y distances between the nodal point of the
% camera in the sceneWorld plane and the point at which the ray will strike
% the plane that contains the camera
slope_xZ =(outputRaySceneWorld3D(2,1)-outputRaySceneWorld3D(1,1))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
slope_yZ =(outputRaySceneWorld3D(2,2)-outputRaySceneWorld3D(1,2))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
cameraPlaneX = outputRaySceneWorld3D(1,1)+((extrinsicTranslationVector(3)-outputRaySceneWorld3D(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld3D(1,2)+((extrinsicTranslationVector(3)-outputRaySceneWorld3D(1,3))*slope_yZ);

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
distance = sqrt((extrinsicTranslationVector(1)-cameraPlaneX)^2 + ...
        (extrinsicTranslationVector(2)-cameraPlaneY)^2 );

end % calcCameraNodeDistanceError3D



%% calcVirtualImageRay
function [outputRayEyeWorld] = calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, traceOpticalSystem)
% Returns the unit vector virtual image ray for the initial depth position
%
% Syntax:
%  [outputRayEyeWorld] = calcVirtualImageRay(eyeWorldPoint, theta_p1p2, theta_p1p3, traceOpticalSystem)
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

% Assemble the virtual image ray
outputRayEyeWorld_p1p2 = traceOpticalSystem(eyeWorldPoint(1), eyeWorldPoint(2), theta_p1p2);
outputRayEyeWorld_p1p3 = traceOpticalSystem(eyeWorldPoint(1), eyeWorldPoint(3), theta_p1p3);

% Adjust the p1 (optical axis) position of the rays to have their initial
% position at the same p1
slope =(outputRayEyeWorld_p1p2(2,2)-outputRayEyeWorld_p1p2(1,2))/(outputRayEyeWorld_p1p2(2,1)-outputRayEyeWorld_p1p2(1,1));
zOffset=outputRayEyeWorld_p1p2(1,1)-eyeWorldPoint(1);
outputRayEyeWorld_p1p2(:,1)=outputRayEyeWorld_p1p2(:,1)-zOffset;
outputRayEyeWorld_p1p2(:,2)=outputRayEyeWorld_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld_p1p3(2,2)-outputRayEyeWorld_p1p3(1,2))/(outputRayEyeWorld_p1p3(2,1)-outputRayEyeWorld_p1p3(1,1));
zOffset=outputRayEyeWorld_p1p3(1,1)-eyeWorldPoint(1);
outputRayEyeWorld_p1p3(:,1)=outputRayEyeWorld_p1p3(:,1)-zOffset;
outputRayEyeWorld_p1p3(:,2)=outputRayEyeWorld_p1p3(:,2)-(zOffset*slope);

outputRayEyeWorld = zeros(2,3);

% Combine the two dimensions into a single, 3D ray
outputRayEyeWorld(1,:) = [outputRayEyeWorld_p1p2(1,1) outputRayEyeWorld_p1p2(1,2) outputRayEyeWorld_p1p3(1,2)];
outputRayEyeWorld(2,:) = [outputRayEyeWorld_p1p2(2,1) outputRayEyeWorld_p1p2(2,2) outputRayEyeWorld_p1p3(2,2)];


end % calcVirtualImageRay



