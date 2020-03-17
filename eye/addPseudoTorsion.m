function eyePoseAdjusted = addPseudoTorsion(eyePose,primaryPosition)
% Add "pseudo" torsion to the eyePose to follow Listing's Law
%
% Syntax:
%  eyePoseAdjusted = addPseudoTorsion(eyePose,primaryPosition)
%
% Description:
%   Head-fixed eye movements generally obey Listing's Law, which observes
%   that there is no change in the "true" torsion of the eye with a change
%   in position. When eye rotations are implemented with quaternions, this
%   property is automatically achieved. In this model eye system, however,
%   we implement eye rotations as a series of rotations following the
%   "Helmholtz coordinates". (The motivation for using a series of
%   rotations instead of quaternions is that it allows us to implement
%   separate rotation centers for each direction of eye rotation).
%
%   To create eye movements that obey Listing's Law, we therefore need to
%   add "pseudo torsion" to the eyePose. This correction, and its
%   theoretical motivation, is described in:
%
%       Wong, Agnes MF. "Listing's law: clinical significance and
%       implications for neural control." Survey of ophthalmology 49.6
%       (2004): 563-575.
%
%   The calculation of pseudoTorsion is performed relative to the "primary
%   position" of the eye, which is the position from which any other eye
%   pose can be achieved by rotation of the eye about a single axis.
%
% Inputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%   primaryPosition       - A 1x2 vector that provides the values for 
%                           [eyeAzimuth, eyeElevation] that defines the
%                           eyePose for which the eye is in primary
%                           position.
%
% Outputs:
%   eyePoseAdjusted       - A 1x4 vector. The eyePose incorporating
%                           pseudoTorsion, resulting in an eyePose that
%                           obeys Listing's Law.
%
% Examples:
%{
    % Test that we recover the value reported by Wong 2004 in Fig 3
    eyePose = [30  30 0 3];
    primaryPosition = [0 0];
    eyePoseAdjusted = addPseudoTorsion(eyePose,primaryPosition);
    assert(max(abs(eyePoseAdjusted(3) - 7.9)));
%}

% Obtain the horizontal and vertical position of the eye relative to the
% primary position of the eye
H = -(eyePose(1)-primaryPosition(1));
V = -(eyePose(2)-primaryPosition(2));

% Calculate the pseudoTorsion, following eq 1 of Wong 2004. This
% calculation must be performed in radians, so there is some conversion
% here.
pseudoTorsion = rad2deg(-(deg2rad(H)*deg2rad(V))/2);

% Create the adjusted eyePose
eyePoseAdjusted = eyePose;
eyePoseAdjusted(3) = eyePoseAdjusted(3) + pseudoTorsion;

end

