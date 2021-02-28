function [navarroD,focalPoint,errors,rayPath1,rayPath2] = calcNavarroD(eye, desiredAccommodation, fieldOrigin, rayIntersectionHeight, cameraMedium, effectiveInfinity)
% Returns the lens accommodation parameter for a desired near focal point
%
% Syntax:
%  navarroD = calcNavarroD(eye, desiredAccommodation, rayHeight, cameraMedium)
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
%   The purpose of this routine is to determine the navarroD parameter that
%   produces the desired accommodative state of a model eye. By default,
%   the routine creates an emmetropic right eye, although this behavior is
%   modified by providing key-value pairs as varargin, which are then
%   passed to the createSceneGeometry fucnction.
%
%   The routine searches over navarroD parameter values until a pair of
%   rays that arise from the near focal point intersect each other on the
%   surface of retina.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   desiredAccommodation  - Scalar. The desired accommodation state of the
%                           eye in diopters. Valid values range from zero
%                           (unaccommodated) to +10. The value sets the
%                           near point of the eye (measured from the
%                           first principal point), where distance (mm) =
%                           1000 / diopters.
%   fieldOrigin           - 1x2 or 2x1 vector that provides the coordinates
%                           in degrees of visual angle (as defined w.r.t
%                           the incident node) of the location from which
%                           the measurement of accommodation will be made.
%                           For example, the position of the fovea in
%                           visual space might be provided. The default
%                           value is [0 0], and thus on the optical axis.
%   rayIntersectionHeight - Scalar. The divergent rays will arrive at the
%                           corneal apex separated by 2x this value.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   navarroD              - Scalar. The parameter "D" that is used in the
%                           Navarro lens shape equations. See the function:
%                               human.lens
%
% Examples:
%{
    % Check that the errors are within tolerance
    eye = modelEyeParameters();
    [navarroD, ~, errors] = calcNavarroD(eye, 0, [0 0]);
    assert(errors(1)<1e-3)
%}
%{
    % Define a default model eye
    eye = modelEyeParameters();
    % Find the navarroD for accommodation to 10 diopters at the fovea
    desiredAccommodation = 10;
    fieldOrigin = eye.landmarks.fovea.degField(1:2);
    [navarroD, focalPoint, errors,rayPath1, rayPath2] = calcNavarroD(eye, desiredAccommodation, fieldOrigin);
    % Create the model eye with this accommodation
    eye = modelEyeParameters('navarroD',navarroD);
    % Show the eye and the converging rays
    opticalSystem = assembleOpticalSystem(eye,...
        'surfaceSetName','mediumToRetina','cameraMedium','air',...
        'opticalSystemNumRows',[]);
    plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'surfaceAlpha',0.05);
    plotOpticalSystem('newFigure',false,'rayPath',rayPath1);
    plotOpticalSystem('newFigure',false,'rayPath',rayPath2);
%}


% Handle missing inputs
if nargin<2
    error('calcNavarroD:invalidArguments','Too few input arguments')
end

if nargin==2
    fieldOrigin = [0, 0];
    rayIntersectionHeight = 0.5;
    cameraMedium = 'air';
    effectiveInfinity = 2000;
end

if nargin==3
    rayIntersectionHeight = 0.5;
    cameraMedium = 'air';
    effectiveInfinity = 2000;
end

if nargin==4
    cameraMedium = 'air';
    effectiveInfinity = 2000;
end

if nargin==5
    effectiveInfinity = 2000;
end

% Generate the optical system
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium,...
    'opticalSystemNumRows',[]);

% The rayOriginDistance is 1000/desiredAccommodation
rayOriginDistance = 1000/desiredAccommodation;

% There are several possible cases for the rayOrigin which require
% different behavior. We detect these here and then act upon the cases
% below.
if isinf(rayOriginDistance) && all(fieldOrigin==0)
    rayCase = 'infOpticalAxis';
end
if isinf(rayOriginDistance) && ~all(fieldOrigin==0)
    rayCase = 'infOffAxis';
end
if ~isinf(rayOriginDistance) && all(fieldOrigin==0)
    rayCase = 'finiteOpticalAxis';
end
if ~isinf(rayOriginDistance) && ~all(fieldOrigin==0)
    rayCase = 'finiteOffAxis';
end

% Handle the cases
switch rayCase
    case 'infOpticalAxis'
        rayOrigin = [effectiveInfinity;0;0];

    case 'infOffAxis'
        ray = quadric.anglesToRay([0;0;0],fieldOrigin(1),fieldOrigin(2));
        rayOrigin = ray(:,2).*effectiveInfinity;

    case 'finiteOpticalAxis'
        
        % Identify the principal point of the system
        principalPoint = calcPrincipalPoint(opticalSystem, rayOriginDistance, rayIntersectionHeight);
        
        % Adjust the rayOrigin to be w.r.t. the principalPoint
        rayOrigin = [rayOriginDistance;0;0]+principalPoint;

    case 'finiteOffAxis'

        % Identify the principal point of the system
        principalPoint = calcPrincipalPoint(opticalSystem, rayOriginDistance, rayIntersectionHeight);
        
        % Find the location of the rayOrigin that is at the requested
        % fieldOrigin. This point will be at a distance w.r.t. the incident
        % node
        rayPath = calcNodalRayFromField(eye,fieldOrigin,rayOriginDistance,cameraMedium);
        inputRay = quadric.coordsToRay(rayPath(:,1:2));
        
        % We need to adjust the initial position of the inputRay, along its
        % ray path, so that the magnitude of the initial Position is the
        % desired rayOriginDistance from the principal point. We have to
        % scale along the ray path so that the initial position maintains
        % its visual field position, as this is defined w.r.t. the incident
        % node. I have tried to find the analytic solution for this but
        % have failed, so here is a brute force approach.
        adjustedOrigin = @(x) inputRay(:,1) + inputRay(:,2).*x - principalPoint;
        adjustObj = @(x) rayOriginDistance - norm(adjustedOrigin(x));
        rayOrigin = adjustedOrigin(fzero(adjustObj,0));
        
end

% Set up the objective for the search. The objective examines divering rays
% that arise from the position implied by the desired accommodation, and
% obtains the internal focal point for those rays. The squared distance of
% this internal focal point from the retinal surface is the error to be
% minimized.
myObj = @(p) objective(p,eye,rayOrigin,rayIntersectionHeight,cameraMedium);

% Define p0 and bounds
p0 = 5;
lb = -5;
ub = 100;

% Options
options = optimset('fmincon');
options.Display = 'off';

% Perform the search
navarroD = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective at the solution
[distanceFocalPointToRetina,focalPoint,raySeparationAtFocalPoint,rayPath1, rayPath2] = ...
    objective(navarroD,eye,rayOrigin,rayIntersectionHeight,cameraMedium);

% Combine the errors
errors = [distanceFocalPointToRetina, raySeparationAtFocalPoint];

% Detect and warn if no accurate solution is found, which is the case for
% some combinations of model eyes and accommodation states.
if distanceFocalPointToRetina > 1e-3
    warnString = ['Cannot accurately accommodate the eye to ' num2str(desiredAccommodation) ' diopters'];
    warning('calcNavarroD:cannotFocus',warnString);
end

end


%% LOCAL FUNCTIONS

function [fVal,focalPoint,raySeparationAtFocalPoint,rayPath1, rayPath2] = objective(navarroD,eye,rayOrigin,rayIntersectionHeight,cameraMedium)

% Update the eye with the specified navarroD for the lens, and make sure
% that the value for accommodation is set to empty to avoid infinite
% recursion in the search.
eye.meta.navarroD = navarroD;
eye.meta.accommodation = [];
eye.lens = human.lens( eye );

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Obtain the internal focal point
[focalPoint, raySeparationAtFocalPoint, rayPath1, rayPath2] = ...
    calcInternalFocalPoint(opticalSystem,rayOrigin,rayIntersectionHeight);

% Evaluate the quadric function for the retina at the focal point. The
% property of the implicit form of the quadric is that it has a value of
% zero for points on the quadric surface.
funcS = quadric.vecToFunc(eye.retina.S);
fVal = funcS(focalPoint(1),focalPoint(2),focalPoint(3))^2;

end


