function pupil = pupil( eye )

% The pupil is an aperture in the iris, centered on the optical axis
pupil.center = [eye.iris.center(1) 0 0];

% The actual pupil of the eye is elliptical. Further, the eccentricity and
% theta of the entrance pupil ellipse changes with pupil dilation:
%
%   Wyatt, Harry J. "The form of the human pupil." Vision Research
%   35.14 (1995): 2021-2036.
%
% Wyatt reported the average ellipse parameters for the entrance pupil
% (with the visual axis aligned with camera axis) under dim and bright
% light conditions. We calculate the corresponding parameters of the actual
% pupil on the optical axis. We then fit a hyperbolic tangent (sigmoidal)
% function to the the eccentricity of the actual pupil as a function of the
% actual pupil radius. The theta values observed by Wyatt were close to
% vertically orientated in the dark, and horizontally oriented in the
% light. We find that a slight tilt away from vertical for the dilated
% pupil allows our model to fit the Mathur 2013 obliquity component
% perfectly. When the actual pupil eccentricity is below zero, the theta is
% set to zero (horizontal), and above zero value it is set to ~pi/2
% (vertical). In the forward model, we take the absolute value of the
% eccentricity returned by the parameters for the actual pupil
% eccentrivity.
%{
    % Observed entrance pupil diameters reported in Wyatt 1995.
    entranceRadius = [3.09/2 4.93/2];
    % Wyatt reported an eccentricity of the pupil of 0.21 under dark
    % conditions. We find that using that value produces model results that
    % disagree with Malthur 2013. We have adopted an upper value of 0.18
    % instead. We also use the convention of a negative eccentricity for a
    % horizontal major axis and a positive eccentricity for vertical.
    entranceEccen = [-0.12 0.18];
    % Prepare scene geometry and eye pose aligned with visual axis
    sceneGeometry = createSceneGeometry();
    % Fix the actual pupil eccentricity at 0
    sceneGeometry.pupil.eccenFcnString = '@(x) 0';
    % Obtain the pupil area in the image for each entrance radius
    % assuming no ray tracing
    sceneGeometry.refraction = [];
    pupilImage = pupilProjection_fwd([-sceneGeometry.eye.axes.visual.degField(1), -sceneGeometry.eye.axes.visual.degField(2), 0, entranceRadius(1)],sceneGeometry);
    actualArea(1) = pupilImage(3);
    pupilImage = pupilProjection_fwd([-sceneGeometry.eye.axes.visual.degField(1), -sceneGeometry.eye.axes.visual.degField(2), 0, entranceRadius(2)],sceneGeometry);
    actualArea(2) = pupilImage(3);
    % Add the ray tracing function to the sceneGeometry
    sceneGeometry = createSceneGeometry();
    % Fix the actual pupil eccentricity at 0
    sceneGeometry.pupil.eccenFcnString = '@(x) 0';
    % Search across actual pupil radii to find the values that match
    % the observed entrance areas.
    myPupilEllipse = @(radius) pupilProjection_fwd([-sceneGeometry.eye.axes.visual.degField(1), -sceneGeometry.eye.axes.visual.degField(2), 0, radius],sceneGeometry);
    myArea = @(ellipseParams) ellipseParams(3);
    myObj = @(radius) (myArea(myPupilEllipse(radius))-actualArea(1)).^2;
    actualRadius(1) = fminunc(myObj, entranceRadius(1));
    myObj = @(radius) (myArea(myPupilEllipse(radius))-actualArea(2)).^2;
    actualRadius(2) = fminunc(myObj, entranceRadius(2));
    % Now find the actual pupil eccentricity that produces the
    % observed entrance pupil eccentricity
    sceneGeometry.eye.pupil.thetas=[0 0];
    place = {'eye' 'pupil' 'eccenFcnString'};
    mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
    myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.axes.visual.degField(1), -sceneGeometry.eye.axes.visual.degField(2), 0, actualRadius(1)],mySceneGeom(eccen));
    myEccen = @(ellipseParams) ellipseParams(4);
    myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(1))).^2;
    actualEccen(1) = -fminsearch(myObj, 0.1);
    sceneGeometry.eye.pupil.thetas = [pi/2, pi/2];
    mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
    myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.axes.visual.degField(1), -sceneGeometry.eye.axes.visual.degField(2), 0, actualRadius(2)],mySceneGeom(eccen));
    myEccen = @(ellipseParams) ellipseParams(4);
    myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(2))).^2;
    actualEccen(2) = fminsearch(myObj, 0.2);
    % We then interpolate the observed values, assuming that the
    % observed values are close to asymptote
    actualRadiusInterp = [actualRadius(1)-.5 actualRadius(1) mean(actualRadius) actualRadius(2) actualRadius(2)+.5];
    actualEccenInterp = [actualEccen(1)/0.96 actualEccen(1) mean(actualEccen) actualEccen(2) actualEccen(2)/0.96];
    % Fit a hand-tuned sigmoidal function
    sigFit = @(scaleX, shiftY, scaleY, x) (tanh((x-mean(actualRadius)).*scaleX)+shiftY)*scaleY;
    fitEccen = fit(actualRadiusInterp',actualEccenInterp',sigFit);
    fprintf('pupil.eccenParams = [-%4.3f %4.3f %4.3f %4.3f];\n',mean(actualRadius),fitEccen.scaleX,fitEccen.shiftY,fitEccen.scaleY);
    % Plot the fit
    figure
    plot(actualRadiusInterp,actualEccenInterp,'kx');
    hold on
    plot(0.5:.1:3,fitEccen(0.5:.1:3),'-r');
%}
% Specify the params and equation that defines the actual pupil ellipse.
% This can be invoked as a function using str2func.
pupil.eccenParams = [-1.847 6.431 0.184 0.113];
pupil.eccenFcnString = sprintf('@(x) (tanh((x+%f).*%f)+%f)*%f',pupil.eccenParams(1),pupil.eccenParams(2),pupil.eccenParams(3),pupil.eccenParams(4));

% The theta values of the actual pupil ellipse for eccentricities less
% than, and greater than, zero.
switch eye.meta.eyeLaterality
    case 'Right'
        pupil.thetas = [0  3/7*pi];
    case 'Left'
        pupil.thetas = [0  4/7*pi];
end

end

