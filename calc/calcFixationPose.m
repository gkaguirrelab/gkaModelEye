function [eyePose,navarroD,visualAngles,initialTargetWorld] = calcFixationPose(fixTargetWorld,stopRadius,varargin)
% Returns the oculomotor pose of the eye to fixate a target
%
% Syntax:
%  [eyePose,navarroD,visualAngles] = calcFixationPose(targetWorldCoords,stopRadius,varargin)
%
% Description
%   Given a target location in world coordinates, this routine reports the
%   parameter values that produce an eye that is fixated upon that point.
%   These parameters include the azimuthal and elevational rotation of the
%   eye and the navarroD accomodation parameter.
%
%   If not defined, the radius of the aperture stop is set to provide an
%   entrance pupil diameter of ~3.5 mm, which tends to produce the
%   highest degree of acuity in normal observers.
%
%   The visual angles of the target (w.r.t. the fovea) are also returned.
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
%       Mapp, Alistair P., and Hiroshi Ono. "The rhino-optical phenomenon:
%       Ocular parallax and the visible field beyond the nose." Vision
%       Research 26.7 (1986): 1163-1165.
%   
%
%
% Inputs:
%   fixTargetWorld        - A 3x1 vector that gives the location of the
%                           fixation target in the world coordinate space
%                           in units of mm.
%   stopRadius            - Scalar. The radius of the aperture stop.
%                           Optional.
%   fixTargetDistance     - Scalar that is the Euclidean distance in mm of 
%                           the fixation target from the origin of the
%                           world coordinate frame. Defaults to 1500.
%
% Optional key/value pairs:
%   None, although varargin are passed to createSceneGeometry
%
% Outputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           This is the pose of the eye for which the line
%                           of sight axis intersects the targetWorldCoords.
%   navarroD              - Scalar. The parameter "D" that is used in the
%                           Navarro lens shape equations. See the function:
%                               human.lens
%   visualAngles          - A 1x2 vector. The visual field position of the
%                           target with respect to the fovea.
%   initialTargetWorld    - The fixation point of the eye prior to
%                           rotation.
%
% Examples:
%{
    fixTargetWorld = [200;100;500];
    [eyePose,navarroD,visualAngles,initialTargetWorld] = calcFixationPose(fixTargetWorld);
    sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    fixTargetEye = convertWorldToEyeCoord(fixTargetWorld);
    initialTargetEye = convertWorldToEyeCoord(initialTargetWorld);
    plot3(fixTargetEye(1),fixTargetEye(2),fixTargetEye(3),'*r');
    plot3(initialTargetEye(1),initialTargetEye(2),initialTargetEye(3),'xb');
%}
    
%% Handle nargin
if nargin == 1
    stopRadius = 1.53;
end

% Force the varargin to include skipping the calculation of angular
% magnification effects from lenses
varargin = [varargin,'skipMagCalc',true];

% There is a small error here, as when the eye is rotated to fixate the
% target, the distance from the corneal apex to the target will be slightly
% different.
fixTargetDistance = norm(fixTargetWorld);

% We will need the target in eye coordinates below
fixTargetEye = convertWorldToEyeCoord(fixTargetWorld);

% If we have not been supplied with a navarroD parameter, calculate it
if ~any(strcmp(varargin,'navarroD'))
    
    % Find the navarroD value that provides for accomodation at the target
    % distance
    navarroD = calcAccommodation(1000/fixTargetDistance, varargin{:});
    
    % Add this to the varargin
    varargin = [varargin,'navarroD',navarroD];
    
end

% Create the sceneGeometry
sceneGeometry = createSceneGeometry(varargin{:},'calcLandmarkFovea',true);

% Find the fixation point along the lineOfSight.
[~,~,~,initialTargetWorld] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance);

% Find the position of the fixation target in units of visual angle w.r.t.
% the fovea. First we have to find the retinal location that receives the
% nodal ray from this point
myObj = @(x) quadric.distancePointRay(fixTargetEye',calcNodalRay(sceneGeometry.eye,[x 0]));
x0 = sceneGeometry.eye.landmarks.vertex.geodetic(1:2);
x = fminsearch(myObj,x0);
Gtarget = [x 0];

% Now find the visual angles of this point w.r.t. the fovea
[~, visualAngles ] = calcVisualAngle(sceneGeometry.eye,sceneGeometry.eye.landmarks.fovea.geodetic,Gtarget);

% Flip the sign on the elevation value to fit the convention that a
% positive elevational rotation moves the eye to look at things that are
% higher up.
visualAngles(2) = -visualAngles(2);

% Find the inverse eye rotation that places the target at the point of best
% focus. This search is done in the eye coordinate space
myObj = @(x) norm(convertWorldToEyeCoord(initialTargetWorld) - rotateEyeCoord(convertWorldToEyeCoord(fixTargetWorld), [x(1) x(2) 0 2], sceneGeometry.eye.rotationCenters, 'inverse'));
x = fminsearch(myObj,visualAngles);

% Prepare the variable to return
eyePose = [x(1) x(2) 0 stopRadius];


end


