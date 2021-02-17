function derivedParams = calcDerivedParams(varargin)
% Derive and save inter-dependent parameters of the eyeModel
%
% Syntax:
%  derivedParams = human.calcDerivedParams()
%
% Description:
%   While many of the model eye parameters are fully set by reference to
%   the literature, some are derived by performing measurements using the
%   model itself. These parameters have a complex, inter-dependent
%   relationship with other derived parameters. This routine gathers in one
%   place the computation of all derived parameters, and saves the
%   resulting values to a file. This way, if there are adjustments to the
%   overall model, this routine can be run to update the derived
%   parameters. Absent a change in the code that computers the eye model,
%   however, execution of this routine will not be needed in general use.
%
% Inputs:
%   none
%
% Outputs:
%   derivedParams         - Structure, with fields for each of the derived
%                           parameter values.
%
% Examples:
%{
    % ETTBSkip -- We don't want to remake the params when we are testing
    % Iteratively run the routine to find the stable set of params
    derivedParams = [];
    % Iteratively call the routine, but do not yet save the params to disk
    for ii=1:2
        derivedParams = human.calcDerivedParams('derivedParams',derivedParams,'derivedParamsPath','')
    end
    % One more iteration, and save the resulting params
    human.calcDerivedParams('derivedParams',derivedParams)
%}


%% input parser
p = inputParser;

% Optional
p.addParameter('derivedParams',[],@(x)(isstruct(x) || isempty(x)));
p.addParameter('derivedParamsPath', ...
    replace(mfilename('fullpath'),mfilename(),'derivedParams.mat'), ...
    @(x)(ischar(x) || isempty(x)));
p.addParameter('verbose',true,@islogical);
p.addParameter('showPlots',true,@islogical);

% parse
p.parse(varargin{:})


%% Initialize derivedParams
% There is a boot-strapping problem as the set of parameters are inter-
% dependent. The params are either passed or set to initial values to
% permit the subsequent calculations to proceed. The ability to pass in the
% initial values supports iterative calls to the routine to allow the
% parameter sets to converge.
if isempty(p.Results.derivedParams)
    derivedParams = struct();
    derivedParams.stopEccenParams = [-1.7523 4.7609 0.1800 0.0973];
    derivedParams.defaultRestingNavarroD = 0.842407226562500;
else
    derivedParams = p.Results.derivedParams;
end


%% Resting accommodation
% The accommodative state of the lens of the eye is set by the navarroD
% parameter. We calculate here the navarroD parameter value that places the
% default model eye at resting accommodation, which is typically assumed to
% be 1.5 diopters.

% Alert the user
if p.Results.verbose
    fprintf('Calculating resting accommodation.\n');
end

derivedParams.defaultRestingNavarroD = calcAccommodation(1.5);


%% Aperture stop ellipticity
% Used in human.stop
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
% disagree with Malthur 2013. We have adopted an upper value of 0.175
% instead. I also use the convention of a negative eccentricity for a
% horizontal major axis and a positive eccentricity for vertical.
entranceEccen = [-0.12 0.175];

% Prepare scene geometry including the fovea
sceneGeometry = createSceneGeometry('derivedParams',derivedParams,'calcLandmarkFovea',true);

% Fix the stop eccentricity at 0 and remove refraction
sg = sceneGeometry;
sg.eye.stop.eccenFcnString = '@(x) 0';
sg.refraction = [];

% Obtain the pupil area in the image for each entrance radius
pupilImage = projectModelEye([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, entranceRadius(1)],sg,'nStopPerimPoints',16);
stopArea(1) = pupilImage(3);
pupilImage = projectModelEye([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, entranceRadius(2)],sg,'nStopPerimPoints',16);
stopArea(2) = pupilImage(3);

% Add the ray tracing function to the sceneGeometry
sg = sceneGeometry;

% Fix the stop eccentricity at 0
sg.eye.stop.eccenFcnString = '@(x) 0';

% Search across stop radii to find the values that match the observed
% entrance areas.
myPupilEllipse = @(radius) projectModelEye([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, radius],sg,'nStopPerimPoints',16);
myArea = @(ellipseParams) ellipseParams(3);
myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(1)).^2;
stopRadius(1) = fminunc(myObj, entranceRadius(1));
myObj = @(radius) (myArea(myPupilEllipse(radius))-stopArea(2)).^2;
stopRadius(2) = fminunc(myObj, entranceRadius(2));

% Set some options for the upcoming fminsearch
opts = optimset('fminsearch');
opts = optimset(opts,'Display','off');

% Now find the stop eccentricity that produces the observed entrance
% pupil eccentricity
sg.eye.stop.thetas=[0 0];
place = {'eye' 'stop' 'eccenFcnString'};
mySceneGeom = @(eccen) setfield(sg,place{:},['@(x) ' num2str(eccen)]);
myPupilEllipse = @(eccen) projectModelEye([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, stopRadius(1)],mySceneGeom(eccen),'nStopPerimPoints',16);
myEccen = @(ellipseParams) ellipseParams(4);
myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(1))).^2;
stopEccen(1) = -fminsearch(myObj, 0.1,opts);
sg.eye.stop.thetas = [pi/2, pi/2];
mySceneGeom = @(eccen) setfield(sg,place{:},['@(x) ' num2str(eccen)]);
myPupilEllipse = @(eccen) projectModelEye([-sg.eye.landmarks.fovea.degField(1), -sg.eye.landmarks.fovea.degField(2), 0, stopRadius(2)],mySceneGeom(eccen),'nStopPerimPoints',16);
myEccen = @(ellipseParams) ellipseParams(4);
myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(2))).^2;
stopEccen(2) = fminsearch(myObj, 0.2,opts);

% We then interpolate the observed values, assuming that the observed
% values are close to asymptote
stopRadiusInterp = [stopRadius(1)-.5 stopRadius(1) mean(stopRadius) stopRadius(2) stopRadius(2)+.5];
stopEccenInterp = [stopEccen(1)/0.96 stopEccen(1) mean(stopEccen) stopEccen(2) stopEccen(2)/0.96];

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','curvefit:fit:noStartPoint');

% Fit a hand-tuned sigmoidal function
sigFit = @(scaleX, shiftY, scaleY, x) (tanh((x-mean(stopRadius)).*scaleX)+shiftY)*scaleY;
fitEccen = fit(stopRadiusInterp',stopEccenInterp',sigFit);

% Restore the warning state
warning(warningState);

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
    xlabel('aperture stop radius [mm]');
    ylabel('aperture stop non-loinear ellipticity');
    title('derivedParams.stopEccenParams');
end


%% Save the derivedParams to disk
if ~isempty(p.Results.derivedParamsPath)
    save(p.Results.derivedParamsPath,'derivedParams')
end

end

