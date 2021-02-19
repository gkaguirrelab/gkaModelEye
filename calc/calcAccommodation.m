function navarroD = calcAccommodation(eye, desiredAccommodation, rayHeight, cameraMedium)
% Returns the lens accommodation parameter for a desired near focal point
%
% Syntax:
%  navarroD = calcAccommodation(eye, desiredAccommodation, rayHeight, cameraMedium)
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
%   rayHeight             - Scalar. Distance of the ray origin from the
%                           optical axis.
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
%}


% Handle missing inputs
if nargin<2
    error('calcAccommodation:invalidArguments','Too few input arguments')
end

if nargin==2
    rayHeight = 1;
    cameraMedium = 'air';
end

if nargin==3
    cameraMedium = 'air';
end

% Set up the objective for the search. The objective examines rays that
% arise from the near focal point implied by the desired accommodation, and
% obtains the internal focal point for those rays. The squared distance of
% this internal focal point from the retinal surface is the error to be
% minimized.
myObj = @(x) objective(x,eye,desiredAccommodation,rayHeight,cameraMedium);
options = optimset('fminsearch');
options.Display = 'off';
[navarroD,fVal] = fminsearch(myObj,5,options);

% Detect and warn if no accurate solution is found, which is the case for
% some combinations of model eyes and accommodation states.
if fVal > 1e-6
    warnString = ['Cannot accurately accommodate the eye to ' num2str(accommodationDiopters) ' diopters'];
    warning('calcAccommodation:cannotFocus',warnString);
end

end


%% LOCAL FUNCTIONS

function [fVal,focalPoint,blur] = objective(navarroD,eye,desiredAccommodation,rayHeight,cameraMedium)

% Update the eye with the specified navarroD for the lens
eye.meta.navarroD = navarroD;
eye.lens = human.lens( eye );

% Obtain the optical system for this eye
opticalSystem = assembleOpticalSystem(eye,...
    'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);

% Obtain the internal focal point
[focalPoint, blur] = calcInternalFocalPoint(opticalSystem,1000/desiredAccommodation,rayHeight);

% Takes as input the vector representation of the quadric surface for the
% retina, and a 3D coordinate. The fVal returned is related to the distance
% of the coordinate from the surface, which is given by evaluating the
% implicit function defined by the quadric for the [x y z] values of the
% coordinate. The property of the implicit form of the quadric is that it
% has a value of zero for points on the quadric surface.
funcS = quadric.vecToFunc(eye.retina.S);
fVal = funcS(focalPoint(1),focalPoint(2),focalPoint(3))^2;

end


