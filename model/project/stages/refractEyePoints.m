function [eyePoints, pointLabels, targetIntersectError] = refractEyePoints(eyePoints,pointLabels,sceneGeometry,eyePose,options)
% Subject eyeWorld points to refraction by the cornea and artificial lenses
%
% Syntax:
%  [eyePoints, pointLabels, targetIntersectError] = refractEyePoints(eyePoints,pointLabels,sceneGeometry,eyePose,options)
%
% Description:
%   This step accounts for the effect of corneal and corrective lens
%   refraction upon the appearance of coordinate points from the eye.
%
% Inputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   p                     - Structure. The structure returned by the
%                           parameter parser in the calling function.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%
% Outputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%   targetIntersectError -  A nx1 vector that contains the distance (in mm)
%                           between the pinhole aperture of the camera and
%                           the intersection of a ray on the camera plane.
%                           The value is nan for points not subject to
%                           refraction by the cornea. All values will be
%                           nan if sceneGeometry.refraction is empty.
%


% Extract some values for clarity in the code that follows
pupilRayFunc = options.pupilRayFunc;
replaceReflectedPoints = options.replaceReflectedPoints;
rayTraceErrorThreshold = options.rayTraceErrorThreshold;
borderSearchPrecision = options.borderSearchPrecision;

% Create a nan vector for the ray trace errors
targetIntersectError = nan(size(eyePoints,1),1);

% Identify the eyePoints subject to refraction by the cornea
refractPointsIdx = find(...
    strcmp(pointLabels,'stopPerimeter')+...
    strcmp(pointLabels,'stopPerimeterFront')+...
    strcmp(pointLabels,'stopPerimeterBack')+...
    strcmp(pointLabels,'irisActualPerimeter'));

% Set the refractFlag if we have a valid refraction field
refractFlag = false;
if isfield(sceneGeometry,'refraction') && ~isempty(pupilRayFunc)
    if ~isempty(sceneGeometry.refraction)
        refractFlag = true;
    end
end

% Perform refraction
if refractFlag

    % Assemble the static args for the findPupilRay
    args = {sceneGeometry.cameraPosition.translation, ...
        sceneGeometry.eye.rotationCenters, ...
        sceneGeometry.refraction.stopToMedium.opticalSystem, ...
        sceneGeometry.refraction.mediumToCamera.opticalSystem};
    
    % Pre-allocate variables to hold the results
    virtualPoints = nan(length(refractPointsIdx),3);
    virtualIntersectError = inf(length(refractPointsIdx),1);
    virtualPointLabels = cell(length(refractPointsIdx),1);
    virtualPointLabels(:) = {''};
    
    % Loop through the eyePoints that are to be refracted
    for ii=1:length(refractPointsIdx)
        
        % Get this eye point
        eyePoint=eyePoints(refractPointsIdx(ii),:);
        
        % Perform the computation using the passed function handle.
        [virtualImageRay, ~, intersectError] = ...
            pupilRayFunc(eyePoint, eyePose, args{:});
        virtualPoint = virtualImageRay(1,:);
        
        % Check if the point has encountered total internal reflection or
        % is a bad ray trace
        retainPoint = true;
                
        if isnan(intersectError) || (intersectError > rayTraceErrorThreshold)
            % The eyePoint did not yield a valid image point. We will not
            % retain the point unless we find a replacement.
            retainPoint = false;
            
            % If this eyePoint is on the stop border, search across smaller
            % stop radii to find a replacement point on the boundary of the
            % pupil perimeter that does make it through.
            if replaceReflectedPoints

                % Find the appropriate center target for the search
                switch pointLabels{refractPointsIdx(ii)}
                    case 'stopPerimeter'
                        centerTarget = sceneGeometry.eye.stop.center;
                    case 'stopPerimeterFront'
                        centerTarget = sceneGeometry.eye.stop.center(1)+sceneGeometry.eye.iris.thickness/2;
                    case 'stopPerimeterBack'
                        centerTarget = sceneGeometry.eye.stop.center(1)-sceneGeometry.eye.iris.thickness/2;
                    case 'irisActualPerimeter'
                        centerTarget = sceneGeometry.eye.iris.center;
                end
                
                % Initialize the searchScalar
                searchScalar = 1.0 - borderSearchPrecision;

                % Perform the search
                stillSearching = true;
                while stillSearching

                    % Move the eyePoint along the vector connecting the eye
                    % point to the stop or iris center
                    shiftedEyePoint = eyePoint - (eyePoint - centerTarget).*(1 - searchScalar);

                    % Subject the shifted point to refraction
                    [virtualImageRay, ~, intersectError] = ...
                        pupilRayFunc(shiftedEyePoint, eyePose, args{:});
                    virtualPoint = virtualImageRay(1,:);

                    % Update the search scalar
                    searchScalar = searchScalar - borderSearchPrecision;
                    if searchScalar <= 0
                        stillSearching = false;
                    end
                    
                    % Test if the newly refracted point meets criterion
                    if intersectError < rayTraceErrorThreshold
                        stillSearching = false;
                        retainPoint = true;
                    end
                end
            else
                retainPoint = false;
            end
        end
        
        % If the ray-trace process yielded a point in the image, store the
        % result
        if retainPoint
            
            % Add the refracted point to the set
            virtualPoints(ii,:) = virtualPoint;
            virtualIntersectError(ii) = intersectError;
            
            % Create a label for the virtual image point
            newPointLabel = pointLabels{refractPointsIdx(ii)};
            newPointLabel = strrep(newPointLabel,'stop','pupil');
            newPointLabel = strrep(newPointLabel,'Actual','');
            virtualPointLabels(ii) = {newPointLabel};
            
        end
    end
    
    % Add the virtual image points to the entire set of eye points. First
    % find those points that had a bad ray trace
    badTraceIdx = any(isnan(virtualPoints),2);
    
    % Remove the bad ray trace points from both the original eyePoint set
    % and the virtual image set
    if any(badTraceIdx)
        virtualPoints(badTraceIdx,:) = [];
        virtualPointLabels(badTraceIdx) = [];
        virtualIntersectError(badTraceIdx) = [];
        
        eyePoints(refractPointsIdx(badTraceIdx),:) = [];
        pointLabels(refractPointsIdx(badTraceIdx)) = [];
        targetIntersectError(refractPointsIdx(badTraceIdx)) = [];
    end
    eyePoints = [eyePoints; virtualPoints];
    pointLabels = [pointLabels; virtualPointLabels];
    targetIntersectError = [targetIntersectError; virtualIntersectError];
    
else
    
    % If there is no refraction, then the pupil is simply the stop. Copy
    % these points over to their new names
    for ii=1:length(refractPointsIdx)
        
        % Get this eyeWorld point
        eyePoint=eyePoints(refractPointsIdx(ii),:);
        
        % Add the refracted point to the set
        eyePoints = [eyePoints; eyePoint];
        targetIntersectError = [targetIntersectError; 0];
        
        % Create a label for the virtual image point
        newPointLabel = pointLabels{refractPointsIdx(ii)};
        newPointLabel = strrep(newPointLabel,'stop','pupil');
        newPointLabel = strrep(newPointLabel,'Actual','');
        pointLabels = [pointLabels; {newPointLabel}];
    end
    
end

end % refractEyePoints