function virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry, varargin )
% Function handles to ray tracing equations
%
% Syntax:
%  virtualImageFuncPointer = compileVirtualImageFunc( sceneGeometry )
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
    sceneGeometry.virtualImageFunc = compileVirtualImageFunc( sceneGeometry );
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


%% input parser
p = inputParser;

% Required
p.addRequired('sceneGeometry',@(x)(isstruct(x) || ischar(x)));

% Optional
p.addParameter('functionDirPath',[],@ischar);
p.addParameter('cleanUpCompileDir',true,@islogical);

% parse
p.parse(sceneGeometry, varargin{:})


%% Variable and path setup
% if the sceneGeometry variable is a character vector, then assume it is
% the path to a sceneGeometry file.
if ischar(sceneGeometry)
    load(sceneGeometry)
end

% Create a directory for the compiled functions
if ~isempty(p.Results.functionDirPath)
    compileDir = p.Results.functionDirPath;
    if ~exist(compileDir,'dir')
        mkdir(compileDir)
    end
else
    compileDir = [];
end


%% traceOpticalSystem
% 2D ray tracing through the cornea and any corrective lenses
% 
% Syntax:
%  outputRay = rayTraceFuncs.traceOpticalSystem(h, theta, z)
%
% Description:
%   The sceneGeometry.opticalSystem specifies the spatial arrangement and
%   refractive indices of a model of the optics of the cornea and any
%   corrective lenses as a set of centered spherical surfaces. This routine
%   uses these parameters to construct a ray-tracing function that takes as
%   input coordinates and angle of a ray arising from a point object and
%   returns as output the unit vector of the ray that emerges from the last
%   optical surface.
%
% Inputs:
%   h                     - Scalar, units of mm. This is the height of 
%                           the object point from the optical axis. In the
%                           corneal system, a value of zero corresponds to
%                           the optical axis that passes through the
%                           corneal apex. In eyeWorld coordinate space,
%                           this corresponds to dimension p2 or p3.
%   theta                 - Scalar, units of radians. This is the angle
%                           that the ray initially makes with the optical
%                           axis of the system. A value of zero is a ray
%                           that is parallel to the optical system, and,
%                           for a ray starting from a point behind the
%                           cornea, a theta in the range [0 pi] corresponds
%                           to a ray that diverges from the optical axis 
%                           towards a positive height.
%   z                     - Scalar, units of mm. The starting position of
%                           the ray along the optical axis. In the cornea
%                           system, a value of zero corresponds to the apex
%                           of the corneal surface, and negative values are
%                           behind the cornea. The pupil plane is at
%                           position -3.7 mm or so. In the eyeWorld
%                           coordinate space, this corresponds to dimension
%                           p1.
%
% Outputs:
%   outputRay             - A 2x2 matrix that describes the vector of the
%                           ray emerging from the system at the last 
%                           surface. outputRay(1,:) contains [position,
%                           height] of the starting point of the output
%                           ray, and outputRay(2,:) is the [position,
%                           height] of the point of termination of the unit
%                           vector.
%

% Create the function
rayTraceFuncs.traceOpticalSystem = traceOpticalSystem(sceneGeometry.opticalSystem, compileDir);
% Add saved function files to path
if ~isempty(p.Results.functionDirPath)
    addpath(compileDir,'-end');
end


%% calcCameraNodeDistanceError2D
% 2D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = rayTraceFuncs.calcCameraNodeDistanceError2D.p1p2(cameraTranslationX, cameraTranslationY, cameraTranslationZ, eyeAzimuthRads, eyeElevationRads, eyeTorsionRads, p1, p2, p3, rotationCenterDepth, theta_p1p2)
%  distance = rayTraceFuncs.calcCameraNodeDistanceError2D.p1p3(cameraTranslationX, cameraTranslationY, cameraTranslationZ, eyeAzimuthRads, eyeElevationRads, eyeTorsionRads, p1, p2, p3, rotationCenterDepth, theta_p1p3)
%
% Description:
%   This function returns the distance between the nodal point of the
%   camera and a ray that has exited from the corneal system. This distance
%   is calculated within the sceneWorld coordinates on an X-Y plane that is
%   positioned at the Z location of the camera. The point of intersection
%   of the ray upon the plane is found, and then the Euclidean distance
%   between this impact point and the nodal point of the camera is
%   returned.
%
%   There are separate function handles for the calculation of a theta
%   with respect to the p1p2 axes of the eyeWorld coordinate system and a
%   theta with respect to the p1p3 axes. Because the eye can be rotated,
%   variations in theta in a single 2D coordinate system in the eyeWorld
%   can nonetheless sweep across the X and Y positions in the nodal plane
%   of the camera in the sceneWorld coordinates.
%
%   Practically, this function is used to find a thetas in p1p2 and p1p3 
%   that minimize the distance between the the intersection point of the
%   ray in the camera plane and the nodal point of the camera. At a
%   distance of zero, the ray would enter the pin hole aperture of the
%   camera and thus produce a point on the resulting image.
%
%   We perform searches separately for thetas in the p1p2 and p1p3 planes
%   to minimize distance. This allows us to use fminsearch over a single
%   variable, which is computationally efficient.
%
% Inputs:
%   cameraTranslationX   
%   cameraTranslationY
%   cameraTranslationZ    - Each is a scalar in units of mm. The three
%                           values correspond to the
%                           extrinsicTranslationVector of the sceneGeometry
%                           structure.
%   eyeAzimuthRads
%   eyeElevationRads
%   eyeTorsionRads        - Each is a scalar in units of radians. These
%                           are the first three values of the eyePoses
%                           vector. Note that the eyePoses vector
%                           specifies rotations in degrees. These values
%                           must be converted to radians before being
%                           passed to this function. This is necessary as
%                           the matlab "sind" function does not accept
%                           symbolic variables.
%   p1, p2 p3             - Each is a scalar in units of mm. These are the
%                           three values for an eyeWorld point.
%   rotationCenterDepth   - Scalar in units of mm. This is the first value
%                           in the field sceneGeometry.eye.rotationCenter
%   theta_p1p2 (or p1p3)  - Scalar in units of radians. The angle WRT the
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
%   distance              - Scalar in units of mm. The Euclidean distance
%                           of the intersection point of the ray on the
%                           Z camera plane from the nodal point of the
%                           camera.

% Create the function
[rayTraceFuncs.calcCameraNodeDistanceError2D.p1p2, rayTraceFuncs.calcCameraNodeDistanceError2D.p1p3] = ...
    calcCameraNodeDistanceError2D(rayTraceFuncs.traceOpticalSystem, compileDir);
rayTraceFuncs.calcCameraNodeDistanceError2D.argumentNames = {'eyeWorldPoint','extrinsicTranslationVector','[eyeAzimuthRads, eyeElevationRads, eyeTorsionRads]','aziRotCenter_p1p2','eleRotCenter_p1p3','torRotCenter_p2p3','theta'};
% Add saved function files to path
if ~isempty(p.Results.functionDirPath)
    addpath(compileDir,'-end');
end


%% calcCameraNodeDistanceError3D
% 3D distance of ray intersection on camera plane from camera node
%
% Syntax:
%  distance = rayTraceFuncs.calcCameraNodeDistanceError3D(cameraTranslationX, cameraTranslationY, cameraTranslationZ, eyeAzimuthRads, eyeElevationRads, eyeTorsionRads, p1, p2, p3, rotationCenterDepth, theta_p1p2, theta_p1p3)
%
% Description:
%   This function is similar to the calcCameraNodeDistanceError2D, except that
%   it takes as input theta in both the p1p2 and p1p3 planes. The distance
%   value that is returned is still the Euclidean distance between the
%   intersection point of the output ray on the Z camera plane and the 
%   nodal point of the camera.
%

% Create the function
rayTraceFuncs.calcCameraNodeDistanceError3D = ...
     calcCameraNodeDistanceError3D(rayTraceFuncs.traceOpticalSystem, compileDir);
% Add saved function files to path
if ~isempty(p.Results.functionDirPath)
    addpath(compileDir,'-end');
end


%% calcVirtualImageRay
% Returns the unit vector virtual image ray for the initial depth position
%
% Syntax:
%  outputRay = rayTraceFuncs.calcVirtualImageRay(p1, p2, p3, theta_p1p2, theta_p1p3)
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
%   p1, p2, p3, theta_p1p2, theta_p1p3 - Interpretation as above for
%                           calcCameraNodeDistanceError2D
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

% Create the function
rayTraceFuncs.calcVirtualImageRay = calcVirtualImageRay(rayTraceFuncs.traceOpticalSystem, compileDir);
% Add saved function files to path
if ~isempty(p.Results.functionDirPath)
    addpath(compileDir,'-end');
end


%% virtualImageFunc
% Returns the virtual image coordinates for a point in eyeWorld space
%
% Syntax:
%  [virtualImagePoint, nodalPointIntersectError] = virtualImageFunc(eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters)
%
% Description:
%   Given a point in eyeWorld space, along with the extrinsic translation
%   vector of the camera, the rotational pose of the eye (in degrees), and
%   the rotation centers of the eye, the routine will return the
%   coordinates in eyeWorld space of the virtual image of the supplied
%   point, at a depth equal to the input point depth.
%
%   This function is either compiled and saved as a mex file, or assembled
%   in memory.
%
% Inputs:
%   eyeWorldPoint         - 1x3 vector, with the axes [p1,p2,p3].
%   extrinsicTranslationVector - 3x1 vector, in mm
%   eyeAzimuth            - Scalar, in degrees
%   eyeElevation          - Scalar, in degrees
%   eyeTorsion            - Scalar, in degrees
%   rotationCenters       - Structure, with the fields azi, ele, tor, each
%                           of which contains a 1x3 vector. This is the
%                           field sceneGeometry.eye.rotationCenters.
%
% Outputs:
%   virtualImagePoint     - 1x3 vector, with the axes [p1,p2,p3], with the
%                           p1 value set equal to the p1 value in the input
%                           eyeWorldPoint.
%   nodalPointIntersectError - Scalar. Distance in mm on the plane of the
%                           nodal point of the camera between the nodal
%                           point at the point of intersection of the
%                           virtual ray from the eyeWorld point.
%

% A functionDirPath has been defined, so we will compile the function as a
% mex file and place it on the path
if ~isempty(p.Results.functionDirPath)
    % Define some argument variables so that the compiler can deduce
    % variable types
    args = {[-3 0 0], sceneGeometry.extrinsicTranslationVector, 0, 0, 0, sceneGeometry.eye.rotationCenters};
    % Change to the compile directory
    initialDir = cd(compileDir);
    % Compile the mex file
    codegen -o virtualImageFuncMex virtualImageFuncPreMex -args args
    % Identify the compiled mex file, the suffix of which will vary
    % depending upon the operating system
    fileLocation = dir('virtualImageFuncMex.*');
    % Clean up the compile dir, if requested
    if p.Results.cleanUpCompileDir
        rmdir('codegen', 's');
        delete('calcTraceOpticalSystem.m');
        delete('calcCameraNodeDistanceError2D_p1p2.m');
        delete('calcCameraNodeDistanceError2D_p1p3.m');
        delete('calcCameraNodeDistanceError3D.m');
        delete('calcVirtualImageRay.m');
    end
    % Refresh the path to add the compiled function
    addpath(compileDir,'-end');
    % Change back to the initial directory
    cd(initialDir);
    % Return the path to the function as the output
    virtualImageFuncPointer.handle = @virtualImageFuncMex;
    virtualImageFuncPointer.path = fullfile(fileLocation.folder,fileLocation.name);
    virtualImageFuncPointer.opticalSystem = sceneGeometry.opticalSystem;
    % Save a copy of this variable in the function directory. The saved
    % variable may be used to re-instantiate the function at a later point.
    filePath = fullfile(fileLocation.folder,'virtualImageFuncPointer');
    save(filePath,'virtualImageFuncPointer');
else
    % Create the function in memory only
    virtualImageFuncPointer.handle = @(eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters) virtualImageFuncMatlab( eyeWorldPoint, extrinsicTranslationVector, eyeAzimuth, eyeElevation, eyeTorsion, rotationCenters, rayTraceFuncs);
    virtualImageFuncPointer.path = [];
    virtualImageFuncPointer.opticalSystem = sceneGeometry.opticalSystem;
end


end % compileVirtualImageFunc -- MAIN



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% traceOpticalSystem
function rayTraceFunc = traceOpticalSystem(opticalSystem, compileDir)
% 2D ray tracing through the cornea and any corrective lenses


% Define some symbolic variables
syms z h theta

% Pass the optical system and the symbolic variables to the ray tracer
outputRay2D = rayTraceCenteredSphericalSurfaces([z h], theta, opticalSystem);

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
if isempty(compileDir)
    rayTraceFunc = matlabFunction(outputRay2D);
else
    functionFileName = fullfile(compileDir,'calcTraceOpticalSystem');
    rayTraceFunc = matlabFunction(outputRay2D,'File',functionFileName);    
end

end % traceOpticalSystem


%% calcCameraNodeDistanceError2D
function [p1p2Func, p1p3Func] = calcCameraNodeDistanceError2D(rayTraceFunc, compileDir)
% 2D distance of ray intersection on camera plane from camera node


% Define some symbolic variables
syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% NOTE:  The rayTraceFunc takes input in the order h, theta, z
outputRayEyeWorld2D_p1p2 = rayTraceFunc(p2,theta_p1p2,p1);
outputRayEyeWorld2D_p1p3 = rayTraceFunc(p3,theta_p1p3,p1);

% Add the third, constant dimension for the output rays
outputRayEyeWorld_p1p2=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) p3;...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) p3];
outputRayEyeWorld_p1p3=[outputRayEyeWorld2D_p1p3(1,1) p2 outputRayEyeWorld2D_p1p3(1,2);...
    outputRayEyeWorld2D_p1p3(2,1) p2 outputRayEyeWorld2D_p1p3(2,2)];

% Prepare to rotate the outputRay into the sceneWorld coordinates
syms cameraTranslationX cameraTranslationY cameraTranslationZ 
syms eyeAzimuthRads eyeElevationRads eyeTorsionRads
syms aziRotCent_p1 aziRotCent_p2
rotationCenterAzi = [aziRotCent_p1 aziRotCent_p2 0];
syms eleRotCent_p1 eleRotCent_p3
rotationCenterEle = [eleRotCent_p1 0 eleRotCent_p3];
syms torRotCent_p2 torRotCent_p3
rotationCenterTor = [0 torRotCent_p2 torRotCent_p3];

RotAzi = [cos(eyeAzimuthRads) -sin(eyeAzimuthRads) 0; sin(eyeAzimuthRads) cos(eyeAzimuthRads) 0; 0 0 1];
RotEle = [cos(eyeElevationRads) 0 sin(eyeElevationRads); 0 1 0; -sin(eyeElevationRads) 0 cos(eyeElevationRads)];
RotTor = [1 0 0; 0 cos(eyeTorsionRads) -sin(eyeTorsionRads); 0 sin(eyeTorsionRads) cos(eyeTorsionRads)];

% Copy eyeWorld rays over the HeadWorld variables
outputRayHeadWorld_p1p2 = outputRayEyeWorld_p1p2;
outputRayHeadWorld_p1p3 = outputRayEyeWorld_p1p3;

% For each of the two coordinates in each ray, shift the eyeWorld ray to
% the rotational center of the eye, rotate for this eye pose, then undo the
% centering

%% Torsion
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenterTor(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotTor*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenterTor(dim);
    end
end
% p1p3
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)-rotationCenterTor(dim);
    end
end
outputRayHeadWorld_p1p3 = (RotTor*(outputRayHeadWorld_p1p3)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)+rotationCenterTor(dim);
    end
end


%% Elevation
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenterEle(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotEle*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenterEle(dim);
    end
end
% p1p3
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)-rotationCenterEle(dim);
    end
end
outputRayHeadWorld_p1p3 = (RotEle*(outputRayHeadWorld_p1p3)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)+rotationCenterEle(dim);
    end
end


%% Azimuth
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)-rotationCenterAzi(dim);
    end
end
outputRayHeadWorld_p1p2 = (RotAzi*(outputRayHeadWorld_p1p2)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p2(coord,dim)=outputRayHeadWorld_p1p2(coord,dim)+rotationCenterAzi(dim);
    end
end
% p1p3
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)-rotationCenterAzi(dim);
    end
end
outputRayHeadWorld_p1p3 = (RotAzi*(outputRayHeadWorld_p1p3)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld_p1p3(coord,dim)=outputRayHeadWorld_p1p3(coord,dim)+rotationCenterAzi(dim);
    end
end

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld_p1p2 = outputRayHeadWorld_p1p2(:,[2 3 1]);
outputRaySceneWorld_p1p3 = outputRayHeadWorld_p1p3(:,[2 3 1]);

% Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
% point at which the ray will strike the plane that contains the camera
slope_xZ =(outputRaySceneWorld_p1p2(2,1)-outputRaySceneWorld_p1p2(1,1))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
slope_yZ =(outputRaySceneWorld_p1p2(2,2)-outputRaySceneWorld_p1p2(1,2))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
cameraPlaneX = outputRaySceneWorld_p1p2(1,1)+((cameraTranslationZ-outputRaySceneWorld_p1p2(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p2(1,2)+((cameraTranslationZ-outputRaySceneWorld_p1p2(1,3))*slope_yZ);

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
if isempty(compileDir)
    p1p2Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p2});
else
    functionFileName = fullfile(compileDir,'calcCameraNodeDistanceError2D_p1p2');
    p1p2Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p2},...
        'File',functionFileName);
end

slope_xZ =(outputRaySceneWorld_p1p3(2,1)-outputRaySceneWorld_p1p3(1,1))/(outputRaySceneWorld_p1p3(2,3)-outputRaySceneWorld_p1p3(1,3));
slope_yZ =(outputRaySceneWorld_p1p3(2,2)-outputRaySceneWorld_p1p3(1,2))/(outputRaySceneWorld_p1p3(2,3)-outputRaySceneWorld_p1p3(1,3));
cameraPlaneX = outputRaySceneWorld_p1p3(1,1)+((cameraTranslationZ-outputRaySceneWorld_p1p3(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p3(1,2)+((cameraTranslationZ-outputRaySceneWorld_p1p3(1,3))*slope_yZ);

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
if isempty(compileDir)
    p1p3Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p3});
else
    functionFileName = fullfile(compileDir,'calcCameraNodeDistanceError2D_p1p3');
    p1p3Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p3},...
            'File',functionFileName);
end


end % calcCameraNodeDistanceError2D


%% calcCameraNodeDistanceError3D
function p1p2p3Func = calcCameraNodeDistanceError3D(rayTraceFunc, compileDir)
% 3D distance of ray intersection on camera plane from camera node


% Define some symbolic variables
syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% NOTE: The rayTraceFunc takes input in the order h, theta, z
outputRayEyeWorld2D_p1p2 = rayTraceFunc(p2,theta_p1p2,p1);
outputRayEyeWorld2D_p1p3 = rayTraceFunc(p3,theta_p1p3,p1);

% Create a 3D output ray system. Shift the p1p3 ray to have the same
% initial p1 value as the p1p2 vector
slope =(outputRayEyeWorld2D_p1p2(2,2)-outputRayEyeWorld2D_p1p2(1,2))/(outputRayEyeWorld2D_p1p2(2,1)-outputRayEyeWorld2D_p1p2(1,1));
zOffset=outputRayEyeWorld2D_p1p2(1,1)-p1;
outputRayEyeWorld2D_p1p2(:,1)=outputRayEyeWorld2D_p1p2(:,1)-zOffset;
outputRayEyeWorld2D_p1p2(:,2)=outputRayEyeWorld2D_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld2D_p1p3(2,2)-outputRayEyeWorld2D_p1p3(1,2))/(outputRayEyeWorld2D_p1p3(2,1)-outputRayEyeWorld2D_p1p3(1,1));
zOffset=outputRayEyeWorld2D_p1p3(1,1)-p1;
outputRayEyeWorld2D_p1p3(:,1)=outputRayEyeWorld2D_p1p3(:,1)-zOffset;
outputRayEyeWorld2D_p1p3(:,2)=outputRayEyeWorld2D_p1p3(:,2)-(zOffset*slope);

% Combine into a single, 3D ray
outputRayEyeWorld3D=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) outputRayEyeWorld2D_p1p3(1,2);...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) outputRayEyeWorld2D_p1p3(2,2)];

% prepare to rotate the outputRay into the sceneWorld coordinates
syms cameraTranslationX cameraTranslationY cameraTranslationZ 
syms eyeAzimuthRads eyeElevationRads eyeTorsionRads
syms aziRotCent_p1 aziRotCent_p2
rotationCenterAzi = [aziRotCent_p1 aziRotCent_p2 0];
syms eleRotCent_p1 eleRotCent_p3
rotationCenterEle = [eleRotCent_p1 0 eleRotCent_p3];
syms torRotCent_p2 torRotCent_p3
rotationCenterTor = [0 torRotCent_p2 torRotCent_p3];

RotAzi = [cos(eyeAzimuthRads) -sin(eyeAzimuthRads) 0; sin(eyeAzimuthRads) cos(eyeAzimuthRads) 0; 0 0 1];
RotEle = [cos(eyeElevationRads) 0 sin(eyeElevationRads); 0 1 0; -sin(eyeElevationRads) 0 cos(eyeElevationRads)];
RotTor = [1 0 0; 0 cos(eyeTorsionRads) -sin(eyeTorsionRads); 0 sin(eyeTorsionRads) cos(eyeTorsionRads)];

% Copy over the outputRay from eye to head world
outputRayHeadWorld3D=outputRayEyeWorld3D;

%% Torsion
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenterTor(dim);
    end
end
outputRayHeadWorld3D = (RotTor*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenterTor(dim);
    end
end


%% Elevation
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenterEle(dim);
    end
end
outputRayHeadWorld3D = (RotEle*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenterEle(dim);
    end
end


%% Azimuth
% p1p2
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)-rotationCenterAzi(dim);
    end
end
outputRayHeadWorld3D = (RotAzi*(outputRayHeadWorld3D)')';
for coord = 1:2
    for dim = 1:3
        outputRayHeadWorld3D(coord,dim)=outputRayHeadWorld3D(coord,dim)+rotationCenterAzi(dim);
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
cameraPlaneX = outputRaySceneWorld3D(1,1)+((cameraTranslationZ-outputRaySceneWorld3D(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld3D(1,2)+((cameraTranslationZ-outputRaySceneWorld3D(1,3))*slope_yZ);

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
if isempty(compileDir)
    p1p2p3Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p2, theta_p1p3});
else
    functionFileName = fullfile(compileDir,'calcCameraNodeDistanceError3D');
    p1p2p3Func = matlabFunction(...
        sqrt((cameraTranslationX-cameraPlaneX)^2 + ...
        (cameraTranslationY-cameraPlaneY)^2 ), ...
        'Vars',{[p1 p2 p3],...
            [cameraTranslationX; cameraTranslationY; cameraTranslationZ],...
            [eyeAzimuthRads, eyeElevationRads, eyeTorsionRads],...
            [aziRotCent_p1, aziRotCent_p2],...
            [eleRotCent_p1, eleRotCent_p3],...
            [torRotCent_p2, torRotCent_p3],...
            theta_p1p2, theta_p1p3},...
        'File',functionFileName);
end

end % calcCameraNodeDistanceError3D


%% calcVirtualImageRay
function [virtualImageRayFunc] = calcVirtualImageRay(rayTraceFunc, compileDir)
% Returns the unit vector virtual image ray for the initial depth position


% Define some symbolic variables
syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% Assemble the virtual image ray
outputRayEyeWorld_p1p2 = rayTraceFunc(p2, theta_p1p2, p1);
outputRayEyeWorld_p1p3 = rayTraceFunc(p3, theta_p1p3, p1);

% Adjust the p1 (optical axis) position of the rays to have the
% their initial position at the same p1
slope =(outputRayEyeWorld_p1p2(2,2)-outputRayEyeWorld_p1p2(1,2))/(outputRayEyeWorld_p1p2(2,1)-outputRayEyeWorld_p1p2(1,1));
zOffset=outputRayEyeWorld_p1p2(1,1)-p1;
outputRayEyeWorld_p1p2(:,1)=outputRayEyeWorld_p1p2(:,1)-zOffset;
outputRayEyeWorld_p1p2(:,2)=outputRayEyeWorld_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld_p1p3(2,2)-outputRayEyeWorld_p1p3(1,2))/(outputRayEyeWorld_p1p3(2,1)-outputRayEyeWorld_p1p3(1,1));
zOffset=outputRayEyeWorld_p1p3(1,1)-p1;
outputRayEyeWorld_p1p3(:,1)=outputRayEyeWorld_p1p3(:,1)-zOffset;
outputRayEyeWorld_p1p3(:,2)=outputRayEyeWorld_p1p3(:,2)-(zOffset*slope);

% Combine the two dimensions into a single, 3D ray
outputRayEyeWorld(1,:) = [outputRayEyeWorld_p1p2(1,1) outputRayEyeWorld_p1p2(1,2) outputRayEyeWorld_p1p3(1,2)];
outputRayEyeWorld(2,:) = [outputRayEyeWorld_p1p2(2,1) outputRayEyeWorld_p1p2(2,2) outputRayEyeWorld_p1p3(2,2)];

% Convert the equation with symbolic variables into a function and either
% return a function handle or save the function to a file
if isempty(compileDir)
    virtualImageRayFunc = matlabFunction(outputRayEyeWorld,...
        'Vars',{p1,p2,p3,theta_p1p2,theta_p1p3});
else
    functionFileName = fullfile(compileDir,'calcVirtualImageRay');
    virtualImageRayFunc = matlabFunction(outputRayEyeWorld,...
        'Vars',{p1,p2,p3,theta_p1p2,theta_p1p3},...
        'File',functionFileName);
end

end % calcVirtualImageRay



