function rotationCenters = rotationCenters( eye )
% Returns the rotationCenters sub-field of an eye model structure
%
% Syntax:
%  rotationCenters = human.rotationCenters( eye )
%
% Description:
%   The centers of rotation of the eye are calculated for the different
%   directions of rotation (elevation, azimuth) and adjusted for the axial
%   length of the eye. The values also differ by eye laterality because of
%   the nasal displacement of the rotation center.
%
%   The rotation center of the eye is often treated as a single, fixed
%   point. A typical assumption is that the center of rotation in
%   emmetropes is 13.3 mm behind the corneal apex:
%
%       Gunter K. vonNoorden, MD; Emilio C. Campos "Binocular Vision and
%       Ocular Motility Theory and Management of Strabismus" American
%       Orthoptic Journal 51.1 (2001): 161-162.
%
%   The source of this value in the cited text is not entirely clear. It
%   appears to be some compromise between the observed centers of rotation
%   that are obtained for azimuthal and elevation rotations. Measurements
%   by Fry & Hill in 1962 and 1963 find that the center of rotation is
%   slightly nasal to the optical axis of the eye, and differs for
%   horizontal and vertical rotations:
%
%       Fry, G. A., and W. W. Hill. "The center of rotation of the
%       eye." Optometry and Vision Science 39.11 (1962): 581-595.
%
%       Fry, Glenn A., and W. W. Hill. "The mechanics of elevating the
%       eye." Optometry and Vision Science 40.12 (1963): 707-716.
%
%   This difference in the apparent horizontal and vertical radii of the
%   eye was subsequently confirmed:
%
%       Hayami, Takehito, Kazunori Shidoji, and Katsuya Matsunaga. "An
%       ellipsoidal trajectory model for measuring the line of sight."
%       Vision research 42.19 (2002): 2287-2293.
%
%   Fry & Hill report that the average azimuthal center of rotation was
%   14.8 mm posterior to the corneal apex (14.7 in the emmetropes), and
%   0.79 mm nasal to visual axis; and the elevation center of rotation was
%   12.2 mm posterior to the corneal apex (12.0 in the emmetropes) and 0.33
%   mm superior. These measurements were made relative to the visual axis
%   of the eye. While the current model is in optical axis coordinates, the
%   effect of this difference is very small (less than 1/100th of a
%   millimeter).
%
%   Note that the Fry & Hill measurements supersede the earlier, Park &
%   Park measurements that claimed substantial translation of the eye
%   during rotation:
%
%       Park, Russell Smith, and George E. Park. "The center of ocular
%       rotation in the horizontal plane." American Journal of
%       Physiology--Legacy Content 104.3 (1933): 545-552.
%
%   The Park & Park result was due to their assumption that all "sight
%   lines" (i.e., rotations of the visual axis of the eye) pass through the
%   same point in space. Fry & Hill report that some subjects (2 of 31)
%   show translation of the eye with rotation. Also, there is a small,
%   transient retraction of the eye following a saccade that we do not
%   attempt to model:
%
%       Enright, J. T. "The aftermath of horizontal saccades: saccadic
%       retraction and cyclotorsion." Vision research 26.11 (1986):
%       1807-1814.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   rotationCenters       - Structure.
%


% Assign the emetropic rotation centers, by eye laterality.
switch eye.meta.eyeLaterality
    case 'Right'
        rotationCenters.azi = [-14.7 0.79 0];
    case 'Left'
        rotationCenters.azi = [-14.7 -0.79 0];
end
rotationCenters.ele = [-12.0 0 0.33];
rotationCenters.tor = [0 0 0];

% Spherical ametropia is correlated with the axial length of the eye. I
% assume here that the center of rotation reflects this change in length.
% Fry & Hill found that azimuthal rotation depth increased by 0.167 mm for
% each negative diopter of spherical refraction, and elevation rotation
% depth by 0.15 mm for each negative diopter. Dick and colleagues (Figure
% 6) found that for each mm of increase in axial length, the center of
% rotation increased by 0.5 mm:
%
%   Dick, Graham L., Bryan T. Smith, and Peter L. Spanos. "Axial length and
%   radius of rotation of the eye." Clinical and Experimental Optometry
%   73.2 (1990): 43-50.
%
% Given that in the Atchison data the axial length of the eye increases by
% 0.27 mm for each negative diopter of spherical ametropic error, this
% would imply a lengthening of the radius of eye rotation by 0.14 mm, which
% is in good agreement with the Fry & Hill observation of 0.15 - 0.167 mm
% of increase.
%
% I scale the azi and ele rotation centers by the ratio of the vitreous
% chamber axial and vertical radii relative to the emmetropic size

% Obtain the radii of the retinal surface for the emmetropic eye
emmetropicEye = eye;
emmetropicEye.meta.sphericalAmetropia = 0;
emmetropicEye.meta.axialLength = [];
emmetropicRetina = human.retina(emmetropicEye);
retinaRadiiEmmetrope = quadric.radii(emmetropicRetina.S)';

% Obtain the radii of the current eye
retinaRadii = quadric.radii(eye.retina.S)';

% Scale the rotation centers
rotationCenters.azi = rotationCenters.azi .* (retinaRadii./retinaRadiiEmmetrope);
rotationCenters.ele = rotationCenters.ele .* (retinaRadii./retinaRadiiEmmetrope);
rotationCenters.tor = rotationCenters.tor .* (retinaRadii./retinaRadiiEmmetrope);

% Assign the primary position of the eye, which is used for calculation of
% "pseudo" torsion in creating eye movements that obey Listing's Law
rotationCenters.primaryPosition = [0 0];

end

