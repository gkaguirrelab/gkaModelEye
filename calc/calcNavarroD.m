function [navarroD,internalFocalPoint,errors,opticalSystem,rayPath1,rayPath2] = calcNavarroD(eye,desiredAccommodation,fieldAngularPosition,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
% Returns the lens accommodation parameter for a desired near focal point
%
% Syntax:
%  [navarroD,internalFocalPoint,errors,opticalSystem,rayPath1,rayPath2] = calcNavarroD(eye,desiredAccommodation,fieldAngularPosition,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)
%
% Description
%   The refractive power of the crystaline lens of the model is a function
%   of the parameter "D" from Navarro's equations. Adjustments to this
%   parameter can be used to set the accommodative state of the model eye
%   so that it has a requested near focal point. The accommodative state of
%   the eye is specified in units of diopters, where the reciprocal of this
%   value gives the distance from the principal point of the optical
%   system to the focal point. 
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
%   desiredAccommodation  - Scalar. The desired accommodation state of the
%                           eye in diopters. Valid values range from zero
%                           (unaccommodated) to +10. The value sets the
%                           near point of the eye (measured from the
%                           first principal point), where distance (mm) =
%                           1000 / diopters.
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
%   navarroD              - Scalar. The parameter "D" that is used in the
%                           Navarro lens shape equations. See the function:
%                               human.lens
%
% Examples:
%{
    % Simple example: focus an eye at infinity along the longitudinal axis
    eye = modelEyeParameters();
    navarroD = calcNavarroD(eye)
%}
%{
    % Focus fovea of the eye at upon a point 100 mm distant
    eye = modelEyeParameters('accommodation',0);
    desiredAccommodation = 1000/100;
    fieldAngularPosition = eye.landmarks.fovea.degField(1:2);
    angleReferenceCoord = eye.landmarks.incidentNode.coords';
    [navarroD,internalFocalPoint,errors,opticalSystem,rayPath1,rayPath2] = calcNavarroD(eye,desiredAccommodation,fieldAngularPosition,angleReferenceCoord);
    % Show the eye and the converging rays
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'surfaceAlpha',0.05);
    plotOpticalSystem('newFigure',false,'rayPath',rayPath1);
    plotOpticalSystem('newFigure',false,'rayPath',rayPath2);
%}


arguments
    eye (1,1) {isstruct}
    desiredAccommodation (1,1) {mustBeNumeric} = 0
    fieldAngularPosition (2,1) {mustBeNumeric} = [0; 0]
    angleReferenceCoord (3,1) {mustBeNumeric} = [0; 0; 0]
    rayIntersectionHeight (1,1) {mustBeNumeric} = 0.25
    effectiveInfinity (1,1) {mustBeNumeric} = 1e4
    cameraMedium = 'air'
end

% The rayOriginDistance is 1000/desiredAccommodation
rayOriginDistance = 1000/desiredAccommodation;

% Set up the objective for the search.
myObj = @(p) objective(p,eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium);

% Define p0 and bounds
p0 = 0;
lb = -5;
ub = 100;

% Options
options = optimset('fmincon');
options.Display = 'off';

% Perform the search
navarroD = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective at the solution
[distanceFocalPointToRetina,internalFocalPoint,raySeparationAtFocalPoint,opticalSystem,rayPath1,rayPath2] = ...
    objective(navarroD,eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium);

% Combine the errors
errors = [distanceFocalPointToRetina,raySeparationAtFocalPoint];

% Detect and warn if no accurate solution is found, which is the case for
% some combinations of model eyes and accommodation states.
if distanceFocalPointToRetina > 1e-3
    warnString = ['Cannot accurately accommodate the eye to ' num2str(desiredAccommodation) ' diopters'];
    warning('calcNavarroD:cannotFocus',warnString);
end

end


%% LOCAL FUNCTIONS

function [fVal,internalFocalPoint,raySeparationAtFocalPoint,opticalSystem,rayPath1,rayPath2] = objective(navarroD,eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium)

% Update the eye with the specified navarroD for the lens, and make sure
% that the value for accommodation is set to empty to avoid infinite
% recursion in the search.
eye.meta.navarroD = navarroD;
eye.meta.accommodation = [];
eye.lens = human.lens(eye);

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Obtain the principal point for the optical system
distanceReferenceCoord = calcPrincipalPoint(opticalSystem);

% Obtain the internal focal point
[internalFocalPoint,raySeparationAtFocalPoint,rayPath1,rayPath2] = ...
    calcInternalFocalPoint(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,rayIntersectionHeight,effectiveInfinity,cameraMedium);

% Evaluate the quadric function for the retina at the focal point. The
% property of the implicit form of the quadric is that it has a value of
% zero for points on the quadric surface.
funcS = quadric.vecToFunc(eye.retina.S);
fVal = funcS(internalFocalPoint(1),internalFocalPoint(2),internalFocalPoint(3))^2;

end


