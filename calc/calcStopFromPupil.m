function stopRadius = calcStopFromPupil(eye,pupilRadius,cameraDistance,cameraMedium)
% Returns the stop radius that produces the observed entrance pupil radius
%
% Syntax:
%  stopRadius = calcStopFromPupil(eye,entranceRadius,cameraMedium)
%
% Description
%   We sometimes wish to simulate the an eye with an entrance pupil of a
%   specified size. To do so, we need to determine the size of the aperture
%   stop that would produce the desired pupil size. This calculation is
%   performed for the stop as viewed by a camera positioned on the
%   longitudinal axis of the optical system.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   pupilRadius           - Scalar. The size in mm of the entrance pupil.
%   cameraDistance        - Scalar. The distance (in mm) of the camera from
%                           the origin of the longitudinal axis
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   stopRadius           - Scalar. The size in mm of the entrance pupil.
%
% Examples:
%{
    % Find the stop radius that corresponds to a 3 mm pupil radius
    eye = modelEyeParameters('spectralDomain','vis');
    pupilRadius = 3;
    stopRadius = calcStopFromPupil(eye,pupilRadius);
    fprintf('A pupil radius of %2.2f mm is produced by a stop of %2.2f mm\n',pupilRadius,stopRadius);
%}


arguments
    eye (1,1) {isstruct}
    pupilRadius (1,1) {mustBeNumeric}
    cameraDistance (1,1) {mustBeNumeric} = 1000
    cameraMedium = 'air'
end

% Place the eye in a sceneGeometry
sceneGeometry = createSceneGeometry(...
    'eye',eye,...
    'cameraTranslation',[0;0;cameraDistance],...
    'cameraMedium',cameraMedium);

% An anonymous function to extract the area from a ellipse parameters
myArea = @(ellipseParams) ellipseParams(3);

% Obtain the area of pupil in the image if the stop was not subject to
% the refractive effects of the cornea
sgNoRefract = sceneGeometry;
sgNoRefract.refraction = [];
stopImage = projectModelEye([0, 0, 0, pupilRadius],sgNoRefract);
stopArea = myArea(stopImage);

% Objective to match observed entrance area.
myPupilEllipse = @(radius) projectModelEye([0, 0, 0, radius],sceneGeometry);
myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;

% Options
options = optimset('fminunc');
options.Display = 'off';

% Search
stopRadius = fminunc(myObj, pupilRadius, options);

end


