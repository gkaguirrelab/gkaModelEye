function medialCanthus = medialCanthus( eye )
% Returns the medial canthus landmark sub-field of an eye model structure
%
% Syntax
%  medialCanthus = human.landmarks.medialCanthus( eye )
%
% Description:
%   Calculates the position of the medial canthus. For these default
%   values, we assume that the axis of the head is aligned with the optic
%   axis of the eye and camera.
%
%   We begin with the width of the palpebral fissure. This can vary
%   substantially by sex and race. We assume a width of 27.5 mm, which was
%   the mean size observed for Black men ages 20-39 years old:
%
%       Price, Kristina M., et al. "Eyebrow and eyelid dimensions: an
%       anthropometric analysis of African Americans and Caucasians."
%       Plastic and reconstructive surgery 124.2 (2009): 615-623.
%
%   In primary position, the palpebral fissure has a slightly greater
%   extent on the medial as compared to lateral side. Based upon anatomical
%   forensic measurements, we set this ratio as 44/56:
%
%       Stephan, Carl N., and Paavi L. Davidson. "The placement of the
%       human eyeball and canthi in craniofacial identification." Journal
%       of Forensic Sciences 53.3 (2008): 612-619.
%
%   Therefore, we distribute the 27.5 mm palpebral fissure width so that
%   the lateral canthus is 12 mm lateral to the corneal apex, and the
%   medial canthus is 15.5 mm medial to the corneal apex.
%
%   The depth of the lateral canthus to the corneal apex is set to 10mm
%   (figure 3):
%
%       van den Bosch, Willem A., Ineke Leenders, and Paul Mulder.
%       "Topographic anatomy of the eyelids, and the effects of sex and
%       age." British journal of ophthalmology 83.3 (1999): 347-352.
%
%   I am unable to find a study of the depth of the medial canthus relative
%   to the corneal apex. My best estimate from axial images of the soft
%   tissues of the orbit is that this distance about half of lateral
%   canthus depth. Therefore, I set this value to 5 mm.
%
%	In most people there is a positive "canthal angle", meaning that the
%	medial canthus is slightly lower on the face than the lateral canthus.
%	This angle varies with age (due to sagging of the soft tissues of the
%	lateral canthus), but in young people the principle source of biometric
%	variation is ethnicity. Faces with "asian" features have a larger,
%	positive canthal angle (9°) in comparison to faces with "caucasian" (4°) or
%	"black" (5°) features. Here, we adopt a value of just over 4°, which
%	corresponds to a 2mm elevation difference between the lateral and
%	medial canthus.
%
%       Rhee, Seung Chul, Kyoung-Sik Woo, and Bongsik Kwon. "Biometric
%       study of eyelid shape and dimensions of different races with
%       references to beauty." Aesthetic plastic surgery 36.5 (2012):
%       1236-1245.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   medialCanthus         - Structure with the subfield coords
%
% Examples:
%{
%}


medialCanthus.coords = [-5 15.5 -1];

switch eye.meta.eyeLaterality
    case 'Right'
        % No change needed
    case 'Left'
        medialCanthus.coords(2) = -medialCanthus.coords(2);
    otherwise
        error('eye laterality not defined')
end


end
