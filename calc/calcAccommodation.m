function [navarroD, fVal, path1, path2] = calcAccommodation(accommodationDiopters, varargin)
% Returns the lens accommodation parameter for a desired near focal point
%
% Syntax:
%  navarroD = calcAccommodation(accommodationDiopters)
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
%   parallel rays intersect each other on the surface of retina.
%
% Inputs:
%  'accommodationDiopters' - Scalar that supplies the accommodation state
%                           of the eye. Valid values range from zero
%                           (unaccommodated) to +10. The value sets the
%                           distance from the princpal point of the eye to
%                           the focal point on the right, where diopters =
%                           1000 / distance(mm).
%
% Optional key/value pairs:
%   None, although varargin are passed to createSceneGeometry
%
% Outputs:
%   navarroD              - Scalar. The parameter "D" that is used in the
%                           Navarro lens shape equations. See the function:
%                               human.lens
%
% Examples:
%{
    navarroD = calcAccommodation(1.5)
%}
%{
    % Plot the eye and the rays
    [navarroD, fVal, path1, path2] = calcAccommodation(1.5);
    sceneGeometry = createSceneGeometry('navarroD',navarroD);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    plotOpticalSystem('newFigure',false,'rayPath',path1);
    plotOpticalSystem('newFigure',false,'rayPath',path2);
    xlim([-25 5]);
%}


% Force the varargin to include skipping the calculation of angular
% magnification effects from lenses
varargin = [varargin,'skipMagCalc',true];


%% Anonymous functions for the eye
% Define these here to use in a subsequent nonlinear search. Each function
% takes as input a candidate Navarro D value.

% A sceneGeometry. The varargin are passed here to createSceneGeometry to
% modify the particulars of the eye.
myScene=@(x) createSceneGeometry('navarroD',x,varargin{:});

% The optical system for a model eye
mySystem=@(x) getfield(myScene(x),'refraction','cameraToRetina','opticalSystem');


%% Anonymous functions for the rays
% Create a pair of rays that arise from the optical axis
intersectHeight = 1;

% The behavior here handles the special case of a desired accommodation of
% zero.
if accommodationDiopters==0
    % The rays are fixed at parallel
    myR1 = @(x) quadric.normalizeRay([100,-1;intersectHeight,0;0,0]);
    myR2 = @(x) quadric.normalizeRay([100,-1;-intersectHeight,0;0,0]);
else
    % The principal point of the optical system.
    myPrincipalPoint = @(x) calcPrincipalPoint(mySystem(x));
    
    % Account for the depth of the principal point in calculating the
    % position from which the rays arise, as the coordinate space is w.r.t.
    % the front corneal surface.
    myRayOrigin = @(x) (1000/accommodationDiopters) - sum(myPrincipalPoint(x).*[1;0;0]);
    
    % Calculate the angle with which the rays diverge from the optical axis
    % such that they will intersect the plane at the intersect height
    myAngle = @(x) rad2deg(atan2(intersectHeight,-myRayOrigin(x)));
    
    % Define the two rays
    myR1 = @(x) quadric.normalizeRay(quadric.anglesToRay([myRayOrigin(x);0;0],myAngle(x),0));
    myR2 = @(x) quadric.normalizeRay(quadric.anglesToRay([myRayOrigin(x);0;0],-myAngle(x),0));
end


%% Anonymous functions internal focal point

% The point of intersection of the rays within the eye
myInternalFocalPoint = @(x) quadric.distanceRays(rayTraceQuadrics(myR1(x), mySystem(x)),rayTraceQuadrics(myR2(x), mySystem(x)));

% A function to return the quadric vector for the retinal surface
myRetina = @(x) getfield(myScene(x),'eye','retina','S')';

% The objective function is the distance by which the focal point within
% the eye misses the retinal surface
myObj = @(x) surfaceDistance(myRetina(x),myInternalFocalPoint(x))^2;


%% Perform the search
[navarroD, fVal] = fminsearch(myObj,5);

% Detect and warn if no accurate solution could be found, which is the
% case for some combinations of model eyes and accommodation states.
if fVal > 1e-6
    warnString = ['Cannot accurately accommodate the eye to ' num2str(accommodationDiopters) ' diopters'];
    warning('calcAccommodation:cannotFocus',warnString);
end

% Obtain the ray paths to return
[~,path1] = rayTraceQuadrics(myR1(navarroD), mySystem(navarroD));
[~,path2] = rayTraceQuadrics(myR2(navarroD), mySystem(navarroD));

end


%% LOCAL FUNCTIONS

function fVal = surfaceDistance(myRetina,coord)
funcS = quadric.vecToFunc(myRetina);
fVal = funcS(coord(1),coord(2),coord(3));
end
