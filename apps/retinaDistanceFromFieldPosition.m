% retinaDistanceFromFieldPosition
%
% This script allows one to calculate the geodesic distance on the retina
% from the fovea for the retinal point that corresponds to a particular
% visual field position (expressed relative to fixation).

% Define the horizontal and vertical visual field position in degrees,
% relative to the (foveated) fixation point
fieldPos = [10 0];

% Define some parameters for the model eye
wavelength = 555;
sphericalAmetropia = -5;
eye = modelEyeParameters('spectralDomain',wavelength,'sphericalAmetropia',sphericalAmetropia);

% Obtain the field position of the fovea w.r.t. the longitudinal axis
fieldPosFovea = eye.landmarks.fovea.degField;

% Some landmarks needed for the calculation
rayOriginDistance = 1000/calcAccommodation(eye);
angleReferenceCoord = eye.landmarks.incidentNode.coords;
distanceReferenceCoord = calcPrincipalPoint(eye);

% Perform the calculation
distance = calcRetinalDistanceFromField(eye,fieldPosFovea,fieldPosFovea+fieldPos,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);

% Report the result
fprintf('The retinal location corresponding to visual field position [%d, %d] is %2.2f mm from the fovea\n',fieldPos,distance);
