function [pointBestFocus, visualAxis, lineOfSight] = calcPointBestFocus(sceneGeometry, stopRadius)
% Returns the point of best focus for a model eye defined in sceneGeometry
%
% Syntax:
%  [pointBestFocus, visualAxis, lineOfSight] = calcPointBestFocus(eye, stopRadius)
%
% Description
%   Given a model eye, the routine identifies the point in space at which
%   the eye has its best focus. This point is defined as the location in
%   space at which the line of sight and visual axis of the eye intersect,
%   or pass with minimum distance. This is also referred to as the "near
%   point" of an eye.
%
%   If not defined, the radius of the aperture stop is set to provide an
%   entrance pupil diameter of ~3.5 mm, which tends to produce the
%   highest degree of acuity in normal observers. The fixation target is
%   assumed to 1500 mm unless set.
%
%   The routine requires that the field sceneGeometry.eye.landmarks.fovea
%   be defined.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   stopRadius            - Scalar. The radius of the aperture stop
%
% Outputs:
%   pointBestFocus        - A 3x1 vector that gives the location of the
%                           point of best focus in eye coordinates
%                           (p1, p2, p3).
%
% Examples:
%{
    % Plot the eye and the visual and line of sight axes
    navarroD = calcAccommodation(10);
    sceneGeometry = createSceneGeometry('navarroD',navarroD,'calcLandmarkFovea',true);
    [pointBestFocus, visualAxis, lineOfSight] = calcPointBestFocus(sceneGeometry);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true);
    plotOpticalSystem('newFigure',false,'rayPath',[visualAxis(:,1) pointBestFocus]);
    plotOpticalSystem('newFigure',false,'rayPath',[lineOfSight(:,1) pointBestFocus]);
%}


% Code to determine the stop radius that corresponds to a pupil diameter of
% 3.5 mm. This value is used as it is found to provide peak acuity for
% normal observers.
%{
    entranceRadius = 3.5/2;
    % Prepare scene geometry and eye pose aligned with visual axis
    sceneGeometry = createSceneGeometry();
    % Obtain the pupil area in the image for the entrance radius
    % assuming no ray tracing
    sceneGeometry.refraction = [];
    pupilImage = projectModelEye([0, 0, 0, entranceRadius],sceneGeometry);
    stopArea = pupilImage(3);
    % Add the ray tracing function to the sceneGeometry
    sceneGeometry = createSceneGeometry();
    % Search across stop radii to find the value that matches the observed
    % entrance area.
    myPupilEllipse = @(radius) projectModelEye([0, 0, 0, radius],sceneGeometry);
    myArea = @(ellipseParams) ellipseParams(3);
    myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
    stopRadius = fminunc(myObj, entranceRadius);
    outline = sprintf('A 3.5mm diameter entrance pupil corresponds to a %2.2fmm stop radius\n',stopRadius);
    fprintf(outline);
%}
    
% Parse inputs
if nargin==1
    stopRadius = 1.53;
end

% Check that the sceneGeometry eye has a foveal landmark
if ~isfield(sceneGeometry.eye.landmarks,'fovea')
    error('The sceneGeometry does not have a fovea defined. Set calcLandmarkFovea to true in the call to createSceneGeometry')
end

% Obtain the visual axis of the eye, which is the nodal ray that intersects
% the fovea
visualAxis = calcNodalRay(sceneGeometry.eye,sceneGeometry.eye.landmarks.fovea.geodetic);

% Obtain the line of sight of the eye focused on effective infinity
lineOfSight = calcLineOfSightRay(sceneGeometry,stopRadius,10000);

% The point of best focus is where these two rays have their closest
% approach
pointBestFocus=quadric.distanceRays(visualAxis,lineOfSight);

end
