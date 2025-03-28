function [eyePose,errors] = calcFixationPose(eye,fieldAngularPosition,targetDistance,stopRadius,addPseudoTorsionFlag,cameraMedium)
% Returns the oculomotor pose of the eye to fixate a target
%
% Syntax:
%  eyePose = calcFixationPose(targetWorldCoords,stopRadius,varargin)
%
% Description
%   Given a target in the visual field (in horizontal and vertical degrees
%   w.r.t. the longitudinal axis of the eye when it is aligned with the
%   camera) and the distance of that target in mm from the incident node of
%   the eye, the routine returns the eye pose parameters of the eye
%   required to place the foveal line-of-sight upon that target.
%
%   Notably, the angles of oculomotor rotation needed to bring the line of
%   sight to the target, and the visual angle between the line of sight and
%   the target, are not equal. This is because the eye rotates about points
%   that are different from the optical center of the eye. This
%   circumstance creates ocular parallax, meaning that the rotation of the
%   eye "shifts" the retinal image, with the degree of this shift
%   influenced by the proximity of the object to the eye.
%
%   For more details see:
%
%       Steinman, R. M., W. B. Cushman, and A. J. Martins. "The precision
%       of gaze." Human neurobiology 1 (1982): 97-109.
%
%       Bingham, Geoffrey P. "Optical flow from eye movement with head
%       immobilized:“Ocular occlusion” beyond the nose." Vision Research
%       33.5-6 (1993): 777-789.
%
%       Mapp, Alistair P., and Hiroshi Ono. "The rhino-optical phenomenon:
%       Ocular parallax and the visible field beyond the nose." Vision
%       Research 26.7 (1986): 1163-1165.
%
%
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   fieldAngularPosition  - 1x2 vector that provides the coordinates in
%                           degrees of visual angle of the target
%                           relative to the longitudinal axis of the eye
%                           when it is aligned with the camera.
%   targetDistance        - Scalar. The distance (in mm) of the origin of
%                           the target from the incident node. Assumed to
%                           be 1500 mm if not defined.
%   stopRadius            - Scalar. Radius of the aperture stop, in mm.
%   addPseudoTorsionFlag  - Logical. Defaults to "true" controls if pseudo-
%                           torsion is added to the eyePose to conform to
%                           Listing's Law. The primary position of the eye
%                           influences the correction, and is found in the
%                           field:
%                               eye.rotationCenters.primaryPosition
%                           For more details see:
%                               /project/stages/addPseudoTorsion.m
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           This is the pose of the eye that places the 
%                           point of regard at the fixation target.
%   errors                - 1x1 matrix with the follow error values:
%                             - L2 norm of the mismatch between the
%                               coordinates of the field target and the
%                               location of the point-of-regard of the eye
%                               following the eyePose rotation.
%
% Examples:
%{
    eye = modelEyeParameters();
    fieldAngularPosition = [-5, 10];
    targetDistance = 1500;
    [eyePose,errors] = calcFixationPose(eye,fieldAngularPosition,targetDistance);
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (1,2) {mustBeNumeric} = [0, 0]
    targetDistance (1,1) {mustBeNumeric} = 1500
    stopRadius (1,1) {mustBeNumeric} = 1.53
    addPseudoTorsionFlag (1,1) {islogical} = true
    cameraMedium = 'air'
end


% Obtain the coordinates of the fovea
rayDestination = eye.landmarks.fovea.coords';

% Derive the line-of-sight for the eye for the specified target distance.
lineOfSightRayPath = calcSightRayToRetina(eye,rayDestination,targetDistance,stopRadius,cameraMedium);

% We retain the incident segment of the line-of-sight ray
lineOfSightRay = quadric.coordsToRay(lineOfSightRayPath(:,1:2));

% Define the coordinates of the desiredFixationPoint, which is defined in terms of angular position w.r.t the incident node
referenceCoord = eye.landmarks.incidentNode.coords';
R = quadric.anglesToRay(referenceCoord,fieldAngularPosition(1),fieldAngularPosition(2));
desiredFixationPoint = R(:,1)+R(:,2).*targetDistance;

% Define the objective
myObj = @(p) objective(p,eye,lineOfSightRay,desiredFixationPoint,addPseudoTorsionFlag);

% p0, which is the angles of the field target minus the position of the
% fovea in field coordinates. This should get us pretty close to the
% solution.
p0 = fieldAngularPosition - eye.landmarks.fovea.degField;

% Bounds
lb = [-90,-90];
ub = [ 90, 90];

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search. The values returned in p are the horizontal and vertical eyePose
% needed to bring the point-of-regard of the eye to the desired fixation
% point.
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Call the objective one more time to get the return values
errors = objective(p,eye,lineOfSightRay,desiredFixationPoint,addPseudoTorsionFlag);

% Assemble a full eyePose vector from p and the passed stopRadius
eyePose = [p(1) p(2) 0 stopRadius];

end


%% Local function
function fVal = objective(p,eye,lineOfSightRay,desiredFixationPoint,addPseudoTorsionFlag)

% The variable "p" holds the candidate horizontal and vertical rotations of
% the eye (in Fick coordinates) relative to [0 0], in which the optical
% axis of the eye and the camera are aligned (i.e., relative to the origin
% of the rotational coordinates). We place these values in an eyePose
% vector, with the last two positions holding the torsion of the eye, and
% the radius of the aperture stop (which is unused here).
eyePose = [p(1), p(2), 0, nan];

% If the addPseudoTorsionFlag is set, then a torsional component is added
% so that the eye movement conforms to Listing's Law.
if addPseudoTorsionFlag
    opts.addPseudoTorsion = addPseudoTorsionFlag;
    sg.eye = eye;
    eyePose = addPseudoTorsion(sg,eyePose,opts);
end

% Apply the eye rotation to the lineOfSightRay
newLineOfSightRay = rotateEyeRay(lineOfSightRay', eyePose, eye.rotationCenters)';

% L2 norm of distance between the desired and obtained fixation point
fVal = norm(quadric.distancePointRay(desiredFixationPoint,newLineOfSightRay));


end

