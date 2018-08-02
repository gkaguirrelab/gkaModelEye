function [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem)
% Returns the position and angle of a resultant ray w.r.t. the optical axis
%
% Syntax:
%  [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem)
%
% Description:
%   This routine implements 3D skew ray tracing through generalized quadric
%   surfaces (including ellipsoids and hyperboloids). The surfaces can have
%   arbitrary positions and orientations.
%
%   This approach to optical systems is described in:
%
%       Langenbucher, Achim, et al. "Ray tracing through a schematic eye
%       containing second-order (quadric) surfaces using 4× 4 matrix
%       notation." Ophthalmic and Physiological Optics 26.2 (2006):
%       180-188.
%
% Inputs:
%   inputRay              - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   opticalSystem         - An mx19 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           state of the ray. Each row contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%                           The first row corresponds to the initial
%                           conditions of the ray. Thus, the refractive
%                           index value given in the first row specifies
%                           the index of the medium in which the ray
%                           arises. The other values for the first row are
%                           ignored. The matrix may have rows of all nans.
%                           These are used to define a fixed sized input
%                           variable for compiled code. They are removed
%                           from the matrix and have no effect.
%
% Outputs:
%   outputRay             - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%
% Examples:
%{
    %% Numerical example in 2D
    % Elagha provides a numerical example in section C which is
    % implemented here as an example.
    %   Elagha, Hassan A. "Generalized formulas for ray-tracing and
    %   longitudinal spherical aberration." JOSA A 34.3 (2017): 335-343.

    % Create the optical system of spherical, aligned surfaces
    boundingBox = [-inf inf -inf inf -inf inf];
    clear opticalSystem
    opticalSystem(1,:)=[nan(1,10) nan nan(1,6) nan 1];
    S = quadric.scale(quadric.unitSphere,10);
    S = quadric.translate(S,[22; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) -1 boundingBox 1 1.2];
    S = quadric.scale(quadric.unitSphere,8);
    S = quadric.translate(S,[9; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1 1];
    S = quadric.scale(quadric.unitSphere,12);
    S = quadric.translate(S,[34; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) -1 boundingBox 1 1.5];
    S = quadric.scale(quadric.unitSphere,10);
    S = quadric.translate(S,[20; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1 1.0];

    % Define an initial ray
    p = [0;0;0];
    u = [1;tand(17.309724);0];
    u = u./sqrt(sum(u.^2));
    inputRay = [p, u];

    % Perform the ray trace
    [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem);

    % Obtain the angle of the ray in the p1p2 plane at each surface.
    % Compare with the Elagha numeric results.
    recoveredThetas = atand(diff(rayPath(2,:))./diff(rayPath(1,:)));
    elaghaThetasDeg = [17.309724 9.479589 4.143784 -5.926743];
    assert(max(abs(recoveredThetas - elaghaThetasDeg))<1e-6);
%}
%{
    %% Pupil point through cornea
    sceneGeometry = createSceneGeometry();
    % Define an initial ray
    p = [sceneGeometry.eye.pupil.center(1); 2; 0];
    u = [1;tand(-15);0];
    u = u./sqrt(sum(u.^2));
    inputRay = [p, u];
    % Perform the ray trace
    [outputRay, rayPath] = rayTraceQuadrics(inputRay, sceneGeometry.refraction.pupilToCamera.opticalSystem);
    % Plot the optical system
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.pupilToCamera,...
        'outputRay',outputRay,'rayPath',rayPath, ...
        'addLighting',true);
%}



%% Initialize variables

% Strip the optical system of any rows which are all nans
opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Determine the number of surfaces
nSurfaces = size(opticalSystem,1);

% Define R (the current state of the ray) as the inputRay
R=inputRay;

% Pre-allocate outputRay and rayPath
outputRay = nan(3,2);
rayPath = nan(3,nSurfaces);
rayPath(:,1)=R(:,1);


%% Peform the ray trace
for ii=2:nSurfaces
    
    % Extract components from the optical system row
    S = quadric.vecToMatrix(opticalSystem(ii,1:10));
    side = opticalSystem(ii,11);
    boundingBox = opticalSystem(ii,12:17);
    mustIntersectFlag = opticalSystem(ii,18);
    nRel = opticalSystem(ii-1,19)/opticalSystem(ii,19);

    % Compute the intersection
    X = quadric.intersectRay(S,R,side,boundingBox);
    
    % Exit if we have missed a "must intersect" surface
    if any(isnan(X)) && mustIntersectFlag
        return
    end
    
    % Get the surface normal. Pass empty for the last input variable to
    % skip checking if the X coordinate is on the surface of the quadric
    N = quadric.surfaceNormal(S,X,side,[]);

    % Get the refracted ray
    R = quadric.refractRay(R,N,nRel);

    % Store the ray path
    rayPath(:,ii) = X;
end


%% Return the output ray
outputRay = R;

end

