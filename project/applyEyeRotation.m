function headPoints = applyEyeRotation(eyePoints,pointLabels,sceneGeometry,p,eyePose)


% Omit the eye rotation centers themselves from rotation.
rotatePointsIdx = find(~contains(pointLabels,{'Rotation'}));

% Copy the eyePoints to the headPoints
headPoints = eyePoints;

% Loop through the points to be rotated. We pass the rotation matrix to
% avoid having to re-calculate this for the rotation of each point.
R = [];
for pp = 1:length(rotatePointsIdx)
    headPoints(rotatePointsIdx(pp),:) = rotateEyeCoord(...
        eyePoints(rotatePointsIdx(pp),:), ...
        eyePose, ...
        sceneGeometry.eye.rotationCenters, ...
        'forward', ...
        R);
end

% If we are projecting a full eye model, label as hidden those posterior
% segment points that are posterior to the most posterior of the centers of
% rotation of the eye, and thus would not be visible to the camera.
if p.Results.fullEyeModelFlag
    seenIdx = strcmp(pointLabels,'retina') .* (headPoints(:,1) >= min([sceneGeometry.eye.rotationCenters.azi(1) sceneGeometry.eye.rotationCenters.ele(1)]));
    seenIdx = logical(seenIdx + ~strcmp(pointLabels,'retina'));
    pointLabels(~seenIdx) = strcat(pointLabels(~seenIdx),'_hidden');
end

end