function derivedParams = calcDerivedParams(varargin)
% Returns the lens sub-field of an eye model structure
%
% Syntax:
%  human.calcDerivedParameters(varargin)
%
% Description:
%   
%
% Inputs:
%   none
%
% Outputs:
%   lens                  - Structure.
%



%% input parser
p = inputParser;

% Optional
p.addParameter('derivedParams',[],@(x)(isstruct(x) || isempty(x)));
p.addParameter('derivedParamsPath',mfilename('fullpath'),@(x)(ischar(x) || isempty(x)));
p.addParameter('verbose',true,@islogical);
p.addParameter('showPlots',true,@islogical);

% parse
p.parse(varargin{:})


%% Initialize derivedParams
% There is a boot-strapping problem as the set of parameters are inter-
% dependent. The params are either passed or set to initial values to
% permit the subsequent calculations to proceed. The ability to pass the
% initial values in supports iterative calls to the routine to allow the
% parameter sets to converge.
if isempty(p.Results.derivedParams)
    derivedParams = struct();
    derivedParams.cornealRotation = [0    2.1347    3.6019];
    derivedParams.accommodationPolyCoef = [0.1448    4.0022   -0.1750];
    derivedParams.stopEccenParams = [-1.752 4.758 0.221 0.102];
else
    derivedParams = p.Results.derivedParams;
end


%% Corneal ellipsoid rotation
% Navarro 2006 reports the rotation of the corneal ellipsoid with respect
% to the keratometric axis, which is the line that connects the fixation
% point of the keratometer (aligned with the instrument optical axis) with
% the center of curvature of the cornea. Here I calculate the rotation of
% the cornea with respect to the optical axis of the eye. The calculation
% is made for an emmetropic eye and assumes that the fixation point of the
% keratometric instrument is at 500 mm from the corneal surface.
% Figure 3 in the Navarro paper reports the sign of the horizontal rotation
% (beta) as being opposite the values that are shown in Table 1 and in the
% text under the results section "optical axis". I adopt the values as
% presented in the figure, as the direction of rotation would otherwise not
% be sensible.

% Alert the user
if p.Results.verbose
    fprintf('Calculating corneal ellipsoid rotation\n');
    tic
end

% Put the fixation target at 500 mm, set the stop radius to a 2mm
% diameter pupil
fixTargetDistance = 500;
stopRadius = 0.8693;
% Obtain the sceneGeometry for an emmetropic eye
sceneGeometry = createSceneGeometry(...
    'sphericalAmetropia',0,...
    'accommodationDiopeters',1000/500,...
    'spectralDomain','vis',...
    'calcLandmarkFovea',true);
% Obtain the fixation angles and fixation target location
[~,~,fixEyePose, fixTargetWorldCoords] = calcLineOfSightRay(sceneGeometry,stopRadius,fixTargetDistance);
% Find the angles of the ray that connects the center of corneal curvature
% with the fixation target
fixTargetEyeCoords = fixTargetWorldCoords([3 1 2]);
cc = sceneGeometry.eye.cornea.front.center;
keratometricAxisRay = quadric.normalizeRay([cc; fixTargetEyeCoords'-cc]');
keratometricAxisAngles = zeros(1,3);
[keratometricAxisAngles(1), keratometricAxisAngles(2)] = quadric.rayToAngles(keratometricAxisRay);
cornealAxisWRTKeratometric = [-2.35 -0.35 0];
% These are the valus OD
cornealAxisWRTOpticalAxis = keratometricAxisAngles + cornealAxisWRTKeratometric;
% Re-arrange to be in terms of the quadric rotations. A rotation that
% directs the corneal apex towards the nasal object space is made about
% the vertical axis.
derivedParams.cornealRotation = fliplr(cornealAxisWRTOpticalAxis);

% Alert the user
if p.Results.verbose
    toc
    fprintf('derivedParams.cornealRotation = [%4.3f %4.3f %4.3f];\n\n', ...
        derivedParams.cornealRotation(1), ...
        derivedParams.cornealRotation(2), ...
        derivedParams.cornealRotation(3));
end




%% Accomodation values
% Because of various imperfections in the model and differences from the
% Navarro paper, it was necessary to "tune" the assigned accomodation
% values so that a requested accomodation state of the emmetropic eye
% results in the expected point of best focus. For example, a non-zero
% setting of the D parameter of the Navarro equations is needed to make the
% emmetropic eye have a point of best focus that approaches infinity.
% I examined the relationship between values of the D parameter and best
% focal distances. After some adjustment by hand, I found that the smallest
% D value that produces the largest best focus distance was ~3.44. Values
% smaller than this cause errors in the model.

% Alert the user
if p.Results.verbose
    fprintf('Calculating Navarro D values for accommodation settings.\n');
    tic
end

navarroD = [3.5, 5, 7.5, 10, 15, 20, 30];
accomodationDiopters = nan(size(navarroD));
sceneGeometry = createSceneGeometry('navarroD',navarroD(1),'calcLandmarkFovea',true);
fovea = sceneGeometry.eye.landmarks.fovea;

for ii = 1:length(navarroD)
    sceneGeometry = createSceneGeometry('navarroD',navarroD(ii));
    sceneGeometry.eye.landmarks.fovea = fovea;
    [outputRayLoS,~] = calcLineOfSightRay(sceneGeometry);
    outputRayVis = calcNodalRay(sceneGeometry.eye,sceneGeometry.eye.landmarks.fovea.geodetic);
    pointOfBestFocus=quadric.distanceRays(outputRayLoS,outputRayVis);
    accomodationDiopters(ii) = 1000/pointOfBestFocus(1);
end

derivedParams.accommodationPolyCoef = polyfit(accomodationDiopters,navarroD,2);

if p.Results.showPlots
    figure
    plot(accomodationDiopters,navarroD,'xk');
    hold on
    plot(1:0.1:10,polyval(derivedParams.accommodationPolyCoef,1:0.1:10),'-r')
end
    
% Alert the user
if p.Results.verbose
    toc
    fprintf('derivedParams.accommodationPolyCoef = [%4.3f %4.3f %4.3f];\n\n', ...
        derivedParams.accommodationPolyCoef(1), ...
        derivedParams.accommodationPolyCoef(2), ...
        derivedParams.accommodationPolyCoef(3));
end




%% Aperture stop ellipticity
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
% (sigmoidal) function to the eccentricity of the stop as a function of the
% stop radius. The theta values observed by Wyatt were close to vertically
% orientated in the dark, and horizontally oriented in the light. We find
% that a slight tilt away from vertical for the dilated pupil allows our
% model to fit the Mathur 2013 obliquity component perfectly. When the stop
% eccentricity is below zero, the theta is set to zero (horizontal), and
% above zero value it is set to ~pi/2 (vertical). In the forward model, we
% take the absolute value of the eccentricity returned by the parameters
% for the stop eccentricity.


% Alert the user
if p.Results.verbose
    fprintf('Calculating aperture stop ellipticity.\n');
    tic
end

% Observed entrance pupil diameters reported in Wyatt 1995.
entranceRadius = [3.09/2 4.93/2];

% Wyatt reported an eccentricity of the pupil of 0.21 under dark
% conditions. I find that using that value produces model results that
% disagree with Malthur 2013. We have adopted an upper value of 0.18
% instead. I also use the convention of a negative eccentricity for a
% horizontal major axis and a positive eccentricity for vertical.
entranceEccen = [-0.12 0.175];

% Prepare scene geometry including the fovea
sceneGeometry = createSceneGeometry('calcLandmarkFovea',true);

% Fix the stop eccentricity at 0 and remove refraction
sg = sceneGeometry;
sg.eye.stop.eccenFcnString = '@(x) 0';
sg.refraction = [];

% Obtain the pupil area in the image for each entrance radius
pupilImage = pupilProjection_fwd([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, entranceRadius(1)],sg,'nStopPerimPoints',16);
stopArea(1) = pupilImage(3);
pupilImage = pupilProjection_fwd([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, entranceRadius(2)],sg,'nStopPerimPoints',16);
stopArea(2) = pupilImage(3);

% Add the ray tracing function to the sceneGeometry
sg = sceneGeometry;

% Fix the stop eccentricity at 0
sg.eye.stop.eccenFcnString = '@(x) 0';

% Search across stop radii to find the values that match the observed
% entrance areas.
myPupilEllipse = @(radius) pupilProjection_fwd([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, radius],sg,'nStopPerimPoints',16);
myArea = @(ellipseParams) ellipseParams(3);
myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
stopRadius(1) = fminunc(myObj, entranceRadius(1));
myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(2)).^2;
stopRadius(2) = fminunc(myObj, entranceRadius(2));

% Now find the stop eccentricity that produces the observed entrance
% pupil eccentricity
sg.eye.stop.thetas=[0 0];
place = {'eye' 'stop' 'eccenFcnString'};
mySceneGeom = @(eccen) setfield(sg,place{:},['@(x) ' num2str(eccen)]);
myPupilEllipse = @(eccen) pupilProjection_fwd([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, stopRadius(1)],mySceneGeom(eccen),'nStopPerimPoints',16);
myEccen = @(ellipseParams) ellipseParams(4);
myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(1))).^2;
stopEccen(1) = -fminsearch(myObj, 0.1);
sg.eye.stop.thetas = [pi/2, pi/2];
mySceneGeom = @(eccen) setfield(sg,place{:},['@(x) ' num2str(eccen)]);
myPupilEllipse = @(eccen) pupilProjection_fwd([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, stopRadius(2)],mySceneGeom(eccen),'nStopPerimPoints',16);
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

% Store the derivedParams
derivedParams.stopEccenParams = [-mean(stopRadius),fitEccen.scaleX,fitEccen.shiftY,fitEccen.scaleY];

if p.Results.verbose
    toc
    fprintf('derivedParams.stopEccenParams = [%4.3f %4.3f %4.3f %4.3f];\n', ...
        derivedParams.stopEccenParams(1), ...
        derivedParams.stopEccenParams(2), ...
        derivedParams.stopEccenParams(3), ...
        derivedParams.stopEccenParams(4));
end

% Plot the fit
if p.Results.showPlots
    figure
    plot(stopRadiusInterp,stopEccenInterp,'kx');
    hold on
    plot(0.5:.1:3,fitEccen(0.5:.1:3),'-r');
end



end

