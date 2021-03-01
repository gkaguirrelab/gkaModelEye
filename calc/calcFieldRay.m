function fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord)
% Ray that arises from the field at a specified angle and distance
%
% Syntax:
%  fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord)
%
% Description
%   This routine is used to define points in object space for various ray
%   tracing operations through the optics of the eye. It is specifically
%   designed to allow specification of a ray origin in terms of angular
%   position with reference to one point (e.g., the incident nodal point)
%   and distance from a different reference point (e.g., the principal
%   point).
%
%   It is useful to define positions in the field in terms of VISUAL ANGLE.
%   This is the horizontal and vertical angle (in degrees) of a point in
%   the field with respect to the approximate incident nodal point of the
%   eye, expressed within the planes of the longitudinal axis. As the
%   approximated nodal characteristic of the eye can change with
%   accommodation, or if measured with ray tracing from different
%   locations, I adopt the standard of measuring the incident nodal point
%   for an eye with the navarroD lens parameter set to zero, and then
%   retaining this as the angleReferenceCoord for all calculations of the
%   visual field. An additional wrinkle is that these angles are expressed
%   with respect to the longitudinal axis of the eye, not the fovea. As
%   the fovea is displaced from the longitudinal axis, the position of the
%   fovea in the visual field is needed to specify visual angle in terms of
%   distance from the point of fixation.
%
%   Sometimes the distance of a point from the eye is important. This is
%   particularly the case when considering the accommodative state of the
%   eye, which is expressed in reciprocal distance (diopters), with the
%   distance expressed relative to the first principal point of the eye.
%   As a consequence, this routine allows one to specify a rayOrigin
%   distance with respect to a point that is different from that used to
%   define angular position.
%
% Inputs:
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in referenceCoord.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the longitudinal axis origin.
%   angleReferenceCoord   - 3x1 vector that provides the coordinate from
%                           which the field angles are calculated. The
%                           incident node is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   distanceReferenceCoord - 3x1 vector that provides the coordinate from
%                           which the rayOriginDistance is calculated. The
%                           The principal point is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%
% Outputs:
%   fieldRay              - 3x2 matrix that specifies the field ray as 
%                           a vector of the form [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%

% First, create a ray that leaves the angleReferenceCoord at the specified
% fieldAngularPositions
Rinitial = quadric.anglesToRay(angleReferenceCoord,fieldAngularPosition(1),fieldAngularPosition(2));

% We now need to determine the location on this ray that is at
% rayOriginDistance from distanceReferenceCoord
myObj = @(d) rayOriginDistance - norm( Rinitial(:,1)+Rinitial(:,2)*d-distanceReferenceCoord );
d = fzero(myObj,rayOriginDistance);

% Assemble the ray with the desired properties
fieldRay = [Rinitial(:,1)+Rinitial(:,2)*d, Rinitial(:,2)];


end