function eyePoseAdjusted = addPseudoTorsion(eyePose,primaryPosition)
% Add "pseudo" torsion to the eyePose to follow Listing's Law
%
% Syntax:
%  eyePoseAdjusted = addPseudoTorsion(eyePose,primaryPosition)
%
% Description:
%   Head-fixed eye movements generally obey Listing's Law, which observes
%   that there is no change in the "true" torsion of the eye (with respect
%   to its optical axis) with a change in position. When eye rotations are
%   implemented with quaternions, this property is automatically achieved.
%   In this model eye system, however, we implement eye rotations as a
%   series of rotations following the "Fick coordinates". (The motivation
%   for using a series of rotations instead of quaternions is that it
%   allows us to implement separate rotation centers for each direction of
%   eye rotation).
%
%   To create eye movements that obey Listing's Law, we therefore need to
%   add "pseudo torsion" to the eyePose. This correction, and its
%   theoretical motivation, is described in:
%
%       Nakayama, K., K. Ciuffreda, and C. Schor. "Kinematics of normal and
%       strabismic eyes." Basic and Clinical Aspects of Binocular Vergence
%       Movements (1983).
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
%}

% Obtain the horizontal and vertical position of the eye relative to the
% primary position of the eye
H = -(eyePose(1)-primaryPosition(1));
V = -(eyePose(2)-primaryPosition(2));

% Calculate the pseudoTorsion, following eq alpha of Nakayama 1983.
A = sind(H)*sind(V);
B = 1+cos(V)*cos(H);

% Need to avoid cases that return complex values
pseudoTorsion = 0;
if B~=0 && abs(A/B)<=1
   pseudoTorsion = -asind(A/B);
end

% Create the adjusted eyePose
eyePoseAdjusted = eyePose;
eyePoseAdjusted(3) = eyePoseAdjusted(3) + pseudoTorsion;

end

