function [eyePoints, pointLabels] = addFullEyeModel(eyePoints, pointLabels, sceneGeometry,p)



%% Define full eye model
% Create a model of the posterior and anterior segments of the eye.

if p.Results.fullEyeModelFlag
    
    % Add points for the stop center, iris center, rotation centers,
    % and retinal landmarks if present
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
    end
    
    % Define points around the perimeter of the iris
    nIrisPerimPoints = p.Results.nIrisPerimPoints;
    perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
    irisPoints = zeros(nIrisPerimPoints,3);
    irisPoints(1:nIrisPerimPoints,3) = ...
        sin(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(3);
    irisPoints(1:nIrisPerimPoints,2) = ...
        cos(perimeterPointAngles)*sceneGeometry.eye.iris.radius + sceneGeometry.eye.iris.center(2);
    irisPoints(1:nIrisPerimPoints,1) = ...
        0 + sceneGeometry.eye.iris.center(1);
    
    % Add the points and labels
    eyePoints = [eyePoints; irisPoints];
    tmpLabels = cell(size(irisPoints,1), 1);
    tmpLabels(:) = {'irisActualPerimeter'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Create the corneal surface points
    corneaPoints = quadric.surfaceGrid(...
        sceneGeometry.eye.cornea.front.S,...
        sceneGeometry.eye.cornea.front.boundingBox,...
        p.Results.corneaMeshDensity, ...
        'parametricPolar');
    
    % Identify the index of the corneal apex
    [~,apexIdx]=max(corneaPoints(:,1));
    
    % Save the corneal apex coordinates
    cornealApex = corneaPoints(apexIdx,:);
    
    % Add the points and labels
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
        p.Results.retinaMeshDensity, ...
        'ellipsoidalPolar');
    
    % Retain those points that are posterior to the iris plane, and have a
    % distance from the optical axis in the p2xp3 plane of greater than the
    % iris radius
    retainIdx = logical(...
        (retinaPoints(:,1) < sceneGeometry.eye.iris.center(1)) .* ...
        sqrt(retinaPoints(:,2).^2+retinaPoints(:,3).^2) > sceneGeometry.eye.iris.radius );
    if all(~retainIdx)
        error('projectModelEye:irisCenterPosition','The iris center is behind the center of the retina');
    end
    retinaPoints = retinaPoints(retainIdx,:);
    
    % Add the points and labels
    eyePoints = [eyePoints; retinaPoints];
    tmpLabels = cell(size(retinaPoints,1), 1);
    tmpLabels(:) = {'retina'};
    pointLabels = [pointLabels; tmpLabels];
    
end

end