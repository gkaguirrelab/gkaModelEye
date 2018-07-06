function [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem)
% Returns the position and angle of a resultant ray w.r.t. the optical axis
%
% Syntax:
%  [outputRay, intersectionCoords, vectorDirections] = rayTraceQuadrics(coordsInitial, inputRay, opticalSystem)
%
% Description:
%   This routine implements 3D skew ray tracing through generalized quadric
%   surfaces (including ellipsoids and hyperboloids). The surfaces can have
%   arbitrary positions and orientations.
%
% Inputs:
%   inputRay              - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   opticalSystem         - An mx18 matrix, where m is the 
%                           number of surfaces in the model, including the
%                           initial state of the ray. Each row contains the
%                           values:
%                               [S side bb n]
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
%                               n     - Refractive index of the surface.
%                           The first row corresponds to the initial
%                           conditions of the ray. Thus, the refractive
%                           index value given in the first row specifies
%                           the index of the medium in which the ray
%                           arises. The center and radius values for the
%                           first row are ignored. The matrix may have rows
%                           of all nans. These are used to define a fixed
%                           sized input variable for compiled code. They
%                           are removed from the matrix and have no effect.
%
% Outputs:
%   outputRay             - 3x2 matrix that specifies the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   rayPath               - An mx3 matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%
% Examples:
%{
    %% Numerical example in 2D
    % This paper provides a numerical example in section C which is
    % implemented here as an example.
    %   Elagha, Hassan A. "Generalized formulas for ray-tracing and
    %   longitudinal spherical aberration." JOSA A 34.3 (2017): 335-343.

    % Create the optical system of spherical, aligned surfaces
    boundingBox = [-inf inf -inf inf -inf inf];
    clear opticalSystem
    opticalSystem(1,:)=[nan(1,10) nan nan(1,6) 1];
    S = quadric.scale(quadric.unitSphere,10);
    S = quadric.translate(S,[22; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) -1 boundingBox 1.2];
    S = quadric.scale(quadric.unitSphere,8);
    S = quadric.translate(S,[9; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1];
    S = quadric.scale(quadric.unitSphere,12);
    S = quadric.translate(S,[34; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) -1 boundingBox 1.5];
    S = quadric.scale(quadric.unitSphere,10);
    S = quadric.translate(S,[20; 0; 0]);
    opticalSystem(end+1,:)=[quadric.matrixToVec(S) 1 boundingBox 1.0];

    % Define an initial ray
    p = [0;0;0];
    u = [1;tand(17.309724);0];
    u = u./sqrt(sum(u.^2));
    inputRay = [p, u];

    % Perform the ray trace
    [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem);

    % Obtain the angle of the ray in the p1p2 plane at each surface.
    % Compare with the Elagha numeric results.
    recoveredThetas = atand(rayPath(:,2,2)./rayPath(:,1,2))';
    elaghaThetasDeg = [17.309724 9.479589 4.143784 -5.926743 -26.583586];
    assert(max(abs(recoveredThetas - elaghaThetasDeg))<1e-4);
%}
%{
    %% Pupil through cornea
    sceneGeometry = createSceneGeometry();
    % Define an initial ray
    p = [sceneGeometry.eye.pupil.center(1); 2; 0];
    u = [1;tand(-15);0];
    u = u./sqrt(sum(u.^2));
    inputRay = [p, u];
    % Perform the ray trace
    outputRay = rayTraceQuadrics(inputRay, sceneGeometry.refraction.opticalSystem);
%}


%% Initialize variables
% Strip the optical system of any rows which are all nans
opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);
% Determine the number of surfaces
nSurfaces = size(opticalSystem,1);
% Define R (the current state of the ray) as the inputRay
R=inputRay;
% Pre-allocate rayPath and normals
rayPath = nan(nSurfaces,3);
rayPath(1,:,:)=R(:,1);


%% Peform the ray trace
for ii=2:nSurfaces
    
    % Extract components from optical system vector
    S = quadric.vecToMatrix(opticalSystem(ii,1:10));
    side = opticalSystem(ii,11);
    boundingBox = opticalSystem(ii,12:17);
    nRel = opticalSystem(ii-1,18)/opticalSystem(ii,18);

    % Compute the intersection, surface normal, and refracted ray
    X = quadric.intersectRay(S,R,side,boundingBox);
    N = quadric.surfaceNormal(S,X,side,[]);
    R = quadric.refractRay(R,N,nRel);

    % Store the ray path and normals
    rayPath(ii,:,:) = X;
end


%% Return the output ray
outputRay = R;

end

