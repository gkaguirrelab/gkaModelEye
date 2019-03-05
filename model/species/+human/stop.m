function stop = stop( eye )
% Returns the stop sub-field of an eye model structure
%
% Syntax:
%  stop = human.stop( eye )
%
% Description:
%   The stop is an aperture in the iris, centered on the optical axis. The
%   stop is modeled as an ellipse, with the eccentricity and theta varying
%   with dilation.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   stop                  - Structure.
%
% 


% Center the stop on the optical axis within the iris
stop.center = [eye.iris.center(1) 0 0];

% The aperture stop of the eye is elliptical. Further, the eccentricity
% and theta of the stop ellipse changes with dilation. The properties of
% the stop are derived from these measurements of the entrance pupil:
%
%   Wyatt, Harry J. "The form of the human stop." Vision Research
%   35.14 (1995): 2021-2036.
%
% Wyatt reported the average ellipse parameters for the entrance pupil
% (with the visual axis aligned with camera axis) under dim and bright
% light conditions. We calculate the corresponding parameters of the
% aperture stop on the optical axis. We then fit a hyperbolic tangent
% (sigmoidal) function to the the eccentricity of the stop as a function of
% the stop radius. The theta values observed by Wyatt were close to
% vertically orientated in the dark, and horizontally oriented in the
% light. We find that a slight tilt away from vertical for the dilated
% pupil allows our model to fit the Mathur 2013 obliquity component
% perfectly. When the stop eccentricity is below zero, the theta is set to
% zero (horizontal), and above zero value it is set to ~pi/2 (vertical). In
% the forward model, we take the absolute value of the eccentricity
% returned by the parameters for the stop eccentricity.
%{
    % Observed entrance pupil diameters reported in Wyatt 1995.
    entranceRadius = [3.09/2 4.93/2];
    % Wyatt reported an eccentricity of the pupil of 0.21 under dark
    % conditions. We find that using that value produces model results that
    % disagree with Malthur 2013. We have adopted an upper value of 0.17
    % instead. We also use the convention of a negative eccentricity for a
    % horizontal major axis and a positive eccentricity for vertical.
    entranceEccen = [-0.12 0.17];
    % Prepare scene geometry and eye pose aligned with line of sight
    sceneGeometry = createSceneGeometry();
    % Fix the stop eccentricity at 0
    sceneGeometry.eye.stop.eccenFcnString = '@(x) 0';
    % Obtain the pupil area in the image for each entrance radius
    % assuming no ray tracing
    sceneGeometry.refraction = [];
    pupilImage = pupilProjection_fwd([-sceneGeometry.eye.axes.lineOfSight.degField(1), -sceneGeometry.eye.axes.lineOfSight.degField(2), 0, entranceRadius(1)],sceneGeometry);
    stopArea(1) = pupilImage(3);
    pupilImage = pupilProjection_fwd([-sceneGeometry.eye.axes.lineOfSight.degField(1), -sceneGeometry.eye.axes.lineOfSight.degField(2), 0, entranceRadius(2)],sceneGeometry);
    stopArea(2) = pupilImage(3);
    % Add the ray tracing function to the sceneGeometry
    sceneGeometry = createSceneGeometry();
    % Fix the stop eccentricity at 0
    sceneGeometry.eye.stop.eccenFcnString = '@(x) 0';
    % Search across stop radii to find the values that match the observed
    % entrance areas.
    myPupilEllipse = @(radius) pupilProjection_fwd([-sceneGeometry.eye.axes.lineOfSight.degField(1), -sceneGeometry.eye.axes.lineOfSight.degField(2), 0, radius],sceneGeometry);
    myArea = @(ellipseParams) ellipseParams(3);
    myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
    stopRadius(1) = fminunc(myObj, entranceRadius(1));
    myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(2)).^2;
    stopRadius(2) = fminunc(myObj, entranceRadius(2));
    % Now find the stop eccentricity that produces the observed entrance
    % pupil eccentricity
    sceneGeometry.eye.stop.thetas=[0 0];
    place = {'eye' 'stop' 'eccenFcnString'};
    mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
    myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.axes.lineOfSight.degField(1), -sceneGeometry.eye.axes.lineOfSight.degField(2), 0, stopRadius(1)],mySceneGeom(eccen));
    myEccen = @(ellipseParams) ellipseParams(4);
    myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(1))).^2;
    stopEccen(1) = -fminsearch(myObj, 0.1);
    sceneGeometry.eye.stop.thetas = [pi/2, pi/2];
    mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
    myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.axes.lineOfSight.degField(1), -sceneGeometry.eye.axes.lineOfSight.degField(2), 0, stopRadius(2)],mySceneGeom(eccen));
    myEccen = @(ellipseParams) ellipseParams(4);
    myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(2))).^2;
    stopEccen(2) = fminsearch(myObj, 0.2);
    % We then interpolate the observed values, assuming that the observed
    % values are close to asymptote
    stopRadiusInterp = [stopRadius(1)-.5 stopRadius(1) mean(stopRadius) stopRadius(2) stopRadius(2)+.5];
    stopEccenInterp = [stopEccen(1)/0.96 stopEccen(1) mean(stopEccen) stopEccen(2) stopEccen(2)/0.96];
    % Fit a hand-tuned sigmoidal function
    sigFit = @(scaleX, shiftY, scaleY, x) (tanh((x-mean(stopRadius)).*scaleX)+shiftY)*scaleY;
    fitEccen = fit(stopRadiusInterp',stopEccenInterp',sigFit);
    fprintf('stop.eccenParams = [-%4.3f %4.3f %4.3f %4.3f];\n',mean(stopRadius),fitEccen.scaleX,fitEccen.shiftY,fitEccen.scaleY);
    % Plot the fit
    figure
    plot(stopRadiusInterp,stopEccenInterp,'kx');
    hold on
    plot(0.5:.1:3,fitEccen(0.5:.1:3),'-r');
%}
% Specify the params and equation that defines the stop ellipse.
% This can be invoked as a function using str2func.
stop.eccenParams = [-1.743 4.784 0.149 0.103];
stop.eccenFcnString = sprintf('@(x) (tanh((x+%f).*%f)+%f)*%f',stop.eccenParams(1),stop.eccenParams(2),stop.eccenParams(3),stop.eccenParams(4));

% The theta values of the stop ellipse for eccentricities less
% than, and greater than, zero.
switch eye.meta.eyeLaterality
    case 'Right'
        stop.thetas = [0  3/7*pi];
    case 'Left'
        stop.thetas = [0  4/7*pi];
end

end

