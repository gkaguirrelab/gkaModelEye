function accommodation = calcAccommodation(eye,fieldAngularPosition,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
% Returns the accommodative state of the eye for a field location
%
% Syntax:
%  accommodation = calcAccommodation(eye,fieldAngularPosition,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
%
% Description
%   The accommodative state of the eye is specified in units of diopters,
%   where the reciprocal of this value gives the distance from the
%   principal point of the optical system to the focal point.
%
%   By default, the calculation is performed with respect to a field
%   position point on the longitudinal axis of the optical system. An
%   alternative (more complicated) choice is to select a
%   fieldAngularPosition and angleReferenceCoord corresponding to the
%   location of the fovea w.r.t. the incidentNode of the eye. A further
%   wrinkle is that the approximation to the incident nodal point shifts
%   with changes in lens properties. Therefore, one convention is to
%   estimate the incident nodal point for an eye accommodated at infinity,
%   and then retain this as the landmark for which visual angle is
%   calculated.

%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in referenceCoord.
%   angleReferenceCoord   - 3x1 vector that provides the coordinate from
%                           which the field angles are calculated. The
%                           incident node is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   rayIntersectionHeight - Scalar. The divergent rays will arrive at the
%                           corneal apex separated by 2x this value.
%   effectiveInfinity     - Scalar. Rays arising from this point or beyond
%                           will be treated as parallel.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   accommodation         - Scalar. The accommodative state (in diopters)
%                           of the eye w.r.t the specified field position.
%
% Examples:
%{
    % Test that we can create and recover a specified accommodation to a
    % point on the longitudinal axis
    desiredAccommodation = 4;
    eye = modelEyeParameters('accommodation',desiredAccommodation);
    measuredAccommodation = calcAccommodation(eye);
    assert(abs(desiredAccommodation-measuredAccommodation)<5e-3)
%}
%{
    % Test that we can create and recover a specified accommodation for the
    % fovea. We first measure the approximation to the incident nodal point
    % in this eye with a navarroD value of zero. We then retain this
    % landmark to define visual field position.
    eye = modelEyeParameters('navarroD',0);
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    fieldAngularPosition = eye.landmarks.fovea.degField;
    desiredAccommodation = 4;
    navarroD = calcNavarroD(eye,desiredAccommodation,fieldAngularPosition,angleReferenceCoord);
    eye = modelEyeParameters('navarroD',navarroD);
    measuredAccommodation = calcAccommodation(eye,fieldAngularPosition,angleReferenceCoord);
    assert(abs(desiredAccommodation-measuredAccommodation)<5e-3)
%}


arguments
    eye (1,1) {isstruct}
    fieldAngularPosition (2,1) {mustBeNumeric} = [0; 0]
    angleReferenceCoord (3,1) {mustBeNumeric} = [0; 0; 0]
    rayIntersectionHeight (1,1) {mustBeNumeric} = 0.25
    effectiveInfinity (1,1) {mustBeNumeric} = 1e4
    cameraMedium = 'air'
end

% Generate the optical system
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium,...
    'opticalSystemNumRows',[]);

% Anonymous function to return the internalFocalPoint based upon the
% rayOriginDistance from the principal point
distanceReferenceCoord = calcPrincipalPoint(opticalSystem);
myFP = @(d) calcInternalFocalPoint(opticalSystem,fieldAngularPosition,d,angleReferenceCoord,distanceReferenceCoord,rayIntersectionHeight,effectiveInfinity);

% The objective is the value of the retinal surface quadric function at the
% focal point
funcS = quadric.vecToFunc(eye.retina.S);
myObj = @(d) objective(myFP(d),funcS);

% Define p0 and bounds
p0 = 1000;
lb = 10;
ub = effectiveInfinity;

% Options
options = optimset('fmincon');
options.Display = 'off';

% Perform the search
rayOriginDistance = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

if rayOriginDistance>=(effectiveInfinity-1e-3)
    accommodation = 0;
else
    accommodation = 1000/rayOriginDistance;
end

end


%% LOCAL FUNCTIONS

% Just need this to distribute the coordinates of the focal point to the
% quadric function
function fVal = objective(internalFocalPoint,funcS)
fVal = funcS(internalFocalPoint(1),internalFocalPoint(2),internalFocalPoint(3))^2;
end


