function [eyePoints, pointLabels] = addFullEyeModel(eyePoints, pointLabels, sceneGeometry,options)
% Define eyeWorld points for the retina, cornea, iris, and other landmarks
%
% Syntax:
%  [eyePoints, pointLabels] = addFullEyeModel(eyePoints, pointLabels, sceneGeometry,p)
%
% Description:
%   Using the contents of the eye field of sceneGeometry, this routine
%   defines a set of coordinate points in eyeWorld space corresponding to
%   quadric surfaces and locations that make up the model eye.
%
% Inputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   p                     - Structure. The structure returned by the
%                           parameter parser in the calling function.
%
% Outputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%


% If this stage is not requested, return
if ~options.fullEyeModelFlag
    return
end


% Add points for the stop center, iris center, rotation centers, and
% landmarks if present
eyePoints = [eyePoints; sceneGeometry.eye.stop.center];
pointLabels = [pointLabels; 'stopCenter'];
eyePoints = [eyePoints; sceneGeometry.eye.iris.center];
pointLabels = [pointLabels; 'irisActualCenter'];
eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.azi];
pointLabels = [pointLabels; 'aziRotationCenter'];
eyePoints = [eyePoints; sceneGeometry.eye.rotationCenters.ele];
pointLabels = [pointLabels; 'eleRotationCenter'];
eyePoints = [eyePoints; 0 0 0];
pointLabels = [pointLabels; 'vertex'];
if isfield(sceneGeometry.eye,'landmarks')
    if isfield(sceneGeometry.eye.landmarks,'fovea')
        eyePoints = [eyePoints; sceneGeometry.eye.landmarks.fovea.coords];
        pointLabels = [pointLabels; 'fovea'];
    end
    if isfield(sceneGeometry.eye.landmarks,'opticDisc')
        eyePoints = [eyePoints; sceneGeometry.eye.landmarks.opticDisc.coords];
        pointLabels = [pointLabels; 'opticDisc'];
    end
    if isfield(sceneGeometry.eye.landmarks,'medialCanthus')
        eyePoints = [eyePoints; sceneGeometry.eye.landmarks.medialCanthus.coords];
        pointLabels = [pointLabels; 'medialCanthus'];
    end
    if isfield(sceneGeometry.eye.landmarks,'lateralCanthus')
        eyePoints = [eyePoints; sceneGeometry.eye.landmarks.lateralCanthus.coords];
        pointLabels = [pointLabels; 'lateralCanthus'];
    end
end

% Define points around the perimeter of the iris
nIrisPerimPoints = options.nIrisPerimPoints;
perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
irisPoints = zeros(nIrisPerimPoints,3);
irisPoints(1:nIrisPerimPoints,3) = ...
    sin(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(3);
irisPoints(1:nIrisPerimPoints,2) = ...
    cos(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(2);
irisPoints(1:nIrisPerimPoints,1) = ...
    0 + sceneGeometry.eye.iris.center(1);

% Add the iris points and labels
eyePoints = [eyePoints; irisPoints];
tmpLabels = cell(size(irisPoints,1), 1);
tmpLabels(:) = {'irisActualPerimeter'};
pointLabels = [pointLabels; tmpLabels];

% Create the corneal surface points
corneaPoints = quadric.surfaceGrid(...
    sceneGeometry.eye.cornea.front.S,...
    sceneGeometry.eye.cornea.front.boundingBox,...
    options.corneaMeshDensity, ...
    'parametricPolar');

% Add the corneal apex, which is the point at ellipsoidal geodetic
% coordinates [0 0 0] (see: quadric.ellipsoidalGeoToCart). The corneal apex
% can be located at either geodetic [0 0 0] or geodetic [0 180 0],
% depending upon the dimensions of the quadric. We test both and take the
% most anterior output.
S = sceneGeometry.eye.cornea.front.S;
cornealApexA = quadric.ellipsoidalGeoToCart([0 0 0],S)';
cornealApexB = quadric.ellipsoidalGeoToCart([0 180 0],S)';
if abs(cornealApexB(1)) < abs(cornealApexA(1))
    cornealApex = cornealApexB;
else
    cornealApex = cornealApexA;
end

% Add the corneal points and labels
eyePoints = [eyePoints; corneaPoints];
tmpLabels = cell(size(corneaPoints,1), 1);
tmpLabels(:) = {'cornea'};
pointLabels = [pointLabels; tmpLabels];

% Add an entry for the corneal apex
eyePoints = [eyePoints; cornealApex];
pointLabels = [pointLabels; 'cornealApex'];

% Create the retina surface vertices
retinaPoints = quadric.surfaceGrid(...
    sceneGeometry.eye.retina.S,...
    sceneGeometry.eye.retina.boundingBox,...
    options.retinaMeshDensity, ...
    'ellipsoidalPolar');

% Retain those points that are posterior to the iris plane, and have a
% distance from the optical axis in the p2xp3 plane of greater than the
% iris radius. This restricts the retinal quadric surface to just those
% points that actually are part of the posterior segment.
retainIdx = logical(...
    (retinaPoints(:,1) < sceneGeometry.eye.iris.center(1)) .* ...
    sqrt(retinaPoints(:,2).^2+retinaPoints(:,3).^2) > sceneGeometry.eye.iris.radius );
if all(~retainIdx)
    error('projectModelEye:irisCenterPosition','The iris center is behind the center of the retina');
end
retinaPoints = retinaPoints(retainIdx,:);

% Add the retina points and labels
eyePoints = [eyePoints; retinaPoints];
tmpLabels = cell(size(retinaPoints,1), 1);
tmpLabels(:) = {'retina'};
pointLabels = [pointLabels; tmpLabels];


end % addFullEyeModel