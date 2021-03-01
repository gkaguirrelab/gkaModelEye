function opticDisc = opticDisc( eye )
% Returns the opticDisc landmark sub-field of an eye model structure
%
% Syntax
%  opticDisc = human.landmarks.opticDisc( eye )
%
% Description:
%   Calculates the position on the retinal surface of the optic disc. This
%   routine requires that the vertex and fovea subfields of landmarks has
%   been defined.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   fovea                 - Structure with the subfields degField,
%                           geodetic, and coords
%



% Check for the vertex subfield
if ~isfield(eye.landmarks,'vertex')
    error('human_landmarks_opticDisc:noVertexDefined','The vertex landmark must be defined in the eye structure');
end

if ~isfield(eye.landmarks,'fovea')
    error('human_landmarks_opticDisc:noFoveaDefined','The fovea landmark must be defined in the eye structure');
end


% The calculation of the position of the retinal landmarks is based upon
% empirical measurements of visual field position. These measurements are
% done in the un-accommodated state, in the visible spectrum, with the eye
% and targets in air. It is also assumed that the position of the retinal
% landmarks were fixed at the point of maturity of the visual system.
% Therefore, the meta values of the passed eye structure are changed here
% to reflect these circumstances for this calculation.
eye.meta.spectralDomain = 'vis';
eye.meta.ageYears = 18;
eye.meta.accommodationDiopters = 0;

% Update the lens field for these values
eye.lens = human.lens(eye);

% The center of the optic disc is nasal and superior to the fovea.
% Rohrschneider measured the distance (in degrees of visual angle) between
% the center of the optic disc and the fovea in 104 people, and found a
% mean horizontal value of 15.5 +- 1.1 degrees SD, and a vertical value of
% 1.5 +- 0.9 degrees SD. There was no dependence of the value upon
% ametropia.
%
%   Rohrschneider, Klaus. "Determination of the location of the fovea on
%   the fundus." Investigative ophthalmology & visual science 45.9 (2004):
%   3257-3258.
%
% By adding the position of fovea in degrees relative to the optical axis,
% we obtain the position of the optic disc in degrees relative to the
% optical axis.
%
switch eye.meta.eyeLaterality
    case 'Right'
        opticDisc.degField = [-15.5 -1.5]+eye.landmarks.fovea.degField;
    case 'Left'
        opticDisc.degField = [15.5 -1.5]+eye.landmarks.fovea.degField;
end

% Now calculate the location on the retina corresponding to this visual
% field location
rayPath = calcNodalRayFromField(eye,opticDisc.degField,1500,eye.landmarks.incidentNode.coords');

% The retinal location is the last point on the rayPath. Store this, and
% obtain the geodetic coordinates
opticDisc.coords = rayPath(:,end)';
opticDisc.geodetic = quadric.cartToEllipsoidalGeo( opticDisc.coords', eye.retina.S )';

end


