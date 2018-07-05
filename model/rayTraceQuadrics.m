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
%   rayPath               - An mx3x2 matrix that provides the refracted ray
%                           at each surface. The value for arrayRay(1,:,:)
%                           is equal to the inputRay, and the value for
%                           arrayRay(end,:,:) is equal to the outputRay.
%
% Examples:
%{
    %% Elagha 2017 numerical example
    % The paper provides a numerical example in section C which is
    % implemented here as an example. Compare the returned theta values
    % with those given on page 340, section C.
    clear figureFlag
    coords = [0 0];
    angleInitial = deg2rad(17.309724);
    figureFlag=true;
    opticalSystem=[nan nan 1; 22 10 1.2; 9 -8 1; 34 12 1.5; 20 -10 1.0];
    [outputRay, angles_p1p2, angles_p1p3, intersectionCoords] = rayTraceEllipsoids(coords, angleInitial, opticalSystem, figureFlag);
    elaghaThetasDeg = [17.309724 9.479589 4.143784 -5.926743 -26.583586];
    assert(max(abs(angles_p1p2' - deg2rad(elaghaThetasDeg)))<1e-5);
%}
%{
    %% Pupil through cornea
    % A model of the passage of a point on the pupil perimeter through
    % the axial cross-section of the cornea (units in mm)
    sceneGeometry = createSceneGeometry();
    [outputRay, angles_p1p2, angles_p1p3, intersectionCoords] = rayTraceEllipsoids([sceneGeometry.eye.pupil.center(1) 2], [deg2rad(-15) 0], sceneGeometry.refraction.opticalSystem, true);
%}
%{
    %% Pupil through cornea, multiple points and rays
    sceneGeometry = createSceneGeometry();
    pupilRadius = 2;
    % Define FigureFlag as a structure, and set the new field to false so
    % that subsequent calls to the ray tracing routine will plot on the
    % same figure. Also, set the textLabels to false to reduce clutter
    figure
    clear figureFlag
    figureFlag.new = false;
    figureFlag.textLabels = false;
    for theta = -25:50:25
        for pupilRadius = -2:4:2
            rayTraceEllipsoids([sceneGeometry.eye.pupil.center(1) pupilRadius], theta, sceneGeometry.refraction.opticalSystem, figureFlag);
        end
    end
%}
%{
    %% Demo non-intersection warning
    coords = [0 0];
    opticalSystem=[nan nan 1.5; 20 10 1.0];
    % This ray will not intersect the surface. The function issues
    % warning and returns an empty outputRay
    theta = deg2rad(45);
    outputRay = rayTraceEllipsoids(coords, theta, opticalSystem);
%}
%{
    %% Demo total internal reflection warning
    coords = [0 0];
    % Make the index of refraction of the surface very high
    opticalSystem=[nan nan 25.0; 20 10 1];
    % This ray encounters total internal reflection. The function issues
    % warning and returns an empty outputRay
    theta = deg2rad(25);
    outputRay = rayTraceEllipsoids(coords, theta, opticalSystem);
%}



%% Initialize variables
% OutputRay set to empty
outputRay = nan(3,2);
% Strip the optical system of any rows which are all nans
opticalSystem=opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);
% Determine the number of surfaces
nSurfaces = size(opticalSystem,1);
% Pre-allocate rayPath
rayPath = nan(nSurfaces,3,2);
rayPath(1,:,:)=inputRay;


%% Peform the ray trace
for ii=2:nSurfaces
    % Extract components from optical system vector
    S = quadric.vecToMatrix(opticalSystem(ii,1:10));
    side = opticalSystem(ii,11);
    boundingBox = opticalSystem(ii,12:17);
    nRel = opticalSystem(ii-1,18)/opticalSystem(ii,18);

    % Compute the intersection, surface normal, and refracted ray
    X = quadric.intersectRay(S,R,side,boundingBox);
    N = quadric.surfaceNormal(S,X,side);
    R = quadric.refractRay(R,N,nRel);

    % Store the ray path
    rayPath(end+1,:,:) = R;
end


end

