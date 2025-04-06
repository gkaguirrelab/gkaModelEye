function imagePoints = projectToImagePlane(worldPoints,sceneGeometry)
% Project the world coordinate points to the image plane
%
% Syntax:
%  imagePoints = projectToImagePlane(worldPoints,sceneGeometry)
%
% Description:
%   Implement the pin-hole camera projection of the set of world 
%   coordinates.
%
%   The image coordinate frame is in units of pixels, and has the
%   dimensions [x, y]:
%
% [0,0]    x
%      +------->
%      |
%   y  |
%      |
%      v
%
% With x being left/right and y being up/down.
%
% Inputs:
%   worldPoints           - nx3 vector. Points in world coordinates.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Outputs:
%   imagePoints           - nx2 vector. Points in image coordinates.
%


% Create the camera position rotation matrix. This is the rotation matrix
% for the position of the camera with respect to the world coordinates. We
% are implementing the Fick rotation system (a.k.a. Tait-Bryan). This
% involves rotation around the axes in the order Z, Y, X. In our
% formulation, this corresponds to applying torsion, elevation, then
% azimuthal rotation
cameraRotationMatrixTors = ...
    [cosd(sceneGeometry.cameraPosition.rotation(3))	-sind(sceneGeometry.cameraPosition.rotation(3))	0; ...
    sind(sceneGeometry.cameraPosition.rotation(3))     cosd(sceneGeometry.cameraPosition.rotation(3))     0; ...
    0                                   0                                       1];

cameraRotationMatrixAzi = ...
    [    1     0     0; ...
    0   cosd(sceneGeometry.cameraPosition.rotation(2))  -sind(sceneGeometry.cameraPosition.rotation(2)); ...
    0   -sind(sceneGeometry.cameraPosition.rotation(2))  cosd(sceneGeometry.cameraPosition.rotation(2))];

cameraRotationMatrixEle = ...
    [cosd(sceneGeometry.cameraPosition.rotation(1))	0 -sind(sceneGeometry.cameraPosition.rotation(1)); ...
    0     1     0; ...
    sind(sceneGeometry.cameraPosition.rotation(1))  0   cosd(sceneGeometry.cameraPosition.rotation(1))];

cameraRotationMatrix = cameraRotationMatrixTors * cameraRotationMatrixEle * cameraRotationMatrixAzi;

% We need to post-multiply the camera rotation matrix by a matrix that
% reverses the Y and Z axes. This is to allow the camera coordinates to
% also be right handed (as are the world coordinates), and to accomodate
% the reversed Y-axis (negative is up) of the image plane.
cameraRotationMatrix = cameraRotationMatrix * [1 0 0; 0 -1 0; 0 0 -1];

% Create the extrinsic matrix. The sceneGeometry structure specifies the
% translation and rotation parameters of camera in world coordinates (i.e.,
% the sceneGeometry structure specifies the camera pose). The extrinsic
% matrix is equal to [R t], where R is transpose of the
% cameraRotationMatrix, and t = -RC, where C is the camerra translation in
% world coordinates
cameraExtrinsicMatrix = ...
    [transpose(cameraRotationMatrix) -transpose(cameraRotationMatrix)*sceneGeometry.cameraPosition.translation];

% The projection matrix is the cameraInstrinsic matrix times the
% cameraExtrinsic matrix.
projectionMatrix = sceneGeometry.cameraIntrinsic.matrix * cameraExtrinsicMatrix;

% What is our total number of points to project?
nEyePoints = size(worldPoints,1);

% Project the world points to camera points, and then to the image plane.
% The world points have a column of ones added to support the
% multiplication with a combined rotation and translation matrix
cameraPoints=(projectionMatrix*[worldPoints, ones(nEyePoints,1)]')';
imagePoints=zeros(nEyePoints,2);
imagePoints(:,1) = ...
    cameraPoints(:,1)./cameraPoints(:,3);
imagePoints(:,2) = ...
    cameraPoints(:,2)./cameraPoints(:,3);

end % projectToImagePlane
