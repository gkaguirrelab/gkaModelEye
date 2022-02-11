function [ fitConeDensityCurcioMm, rawCellDensity ] = cellConeDensityCurcio( dataFileName )
% Brief one line description of the function
%
% Syntax:
%   [ fitTotalConeDensityCurcioMm, rawCellDensity ] = cellConeDensityCurcio( dataFileName )
%
% Description:
%   This routine returns a function variable that provides total RGC
%   density as a function of theta and eccentricity in mm. The fit is
%   derived from the data presented in Curcio et al (1990). Fits to the
%   cardinal meridians are obtained. Interpolation over the parameters of
%   the fit are used to produce a function that returns RGC density for an
%   arbitrary meridian angle.
%
% Inputs:
%   theta            - The desired angle of the density function on
%                           the retinal field. (0=nasal; 90=superior;
%                           180=temporal; 270=inferior)
%
% Optional key/value pairs:
%  'cardinalMeridianAngles' - The polar angles corresponding to the
%                           cardinal medians
%  'splineKnots'          - The number of spline knots to use in the fit
%  'splineOrder'          - The polynomial order of the spline
%  'rgcDensityDataFileName' - The filename of the RGC density file to be
%                           passed to loadRawRGCDensityByEccen. If set to
%                           empty, then the default setting in the load
%                           routine will be used.
%  'makePlots'            - Do we make a figure?
%
% Outputs:
%   fitConeDensityCurcioMm - Anonymous function of theta and
%                           eccentricity (in degrees and mm)
%   rawCellDensity        - Structure that contains the original Curcio
%                           data that was the basis of the fit.
%
% Examples:
%{
    cardinalMeridianAngles = [0, 90, 180, 270];
    cardinalMeridianNames = {'nasal','superior','temporal','inferior'};
    [ fitConeDensityCurcioMm, curcioRawConeDensity ] = cellConeDensityCurcio( );
    fitSupport = 0:0.1:20;
    figure
    for mm = 1:length(cardinalMeridianNames);
        subplot(2,2,mm);
        semilogy(curcioRawConeDensity.support,curcioRawConeDensity.(cardinalMeridianNames{mm}),'xk');
        hold on    
        semilogy(fitSupport,fitConeDensityCurcioMm(cardinalMeridianAngles(mm),fitSupport),'-r');
        xlabel('Eccentricity [mm]'); ylabel('Total cone density per sq mm');
        title(cardinalMeridianNames{mm});
    end
%}
%{
    % Plot the ratio of RGC to cone density across eccentricity
    cardinalMeridianAngles = [0, 90, 180, 270];
    cardinalMeridianNames = {'nasal','superior','temporal','inferior'};
    fitConeDensityCurcioMm = cellConeDensityCurcio( );
    fitTotalRGCDensityCurcioMm = cellTotalRGCDensityCurcio( );
    support = 2:0.1:15;
    figure
    for mm = 1:length(cardinalMeridianNames)
        semilogy(support,fitTotalRGCDensityCurcioMm(cardinalMeridianAngles(mm),support)./fitConeDensityCurcioMm(cardinalMeridianAngles(mm),support),'-');
        hold on
    end
    legend(cardinalMeridianNames);
    xlabel('Eccentricity [mm]'); ylabel('RGC:cone ratio');
%}


arguments
    dataFileName {mustBeTextScalar} = 'curcioRawConeDensity_reportedAverage.mat'
end

% Define the data file to load
fullDataFilePath = fullfile(fileparts(mfilename('fullpath')),'Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography',dataFileName);

% Load and sanity check
tmpLoader = load(fullDataFilePath);
fieldNames = fields(tmpLoader);
if length(fieldNames) ~= 1
    error('Please specify a raw datafile that contains a single structure variable');
end
rawCellDensity = tmpLoader.(fieldNames{1});

cardinalMeridianAngles = [0, 90, 180, 270];
cardinalMeridianNames = {'nasal','superior','temporal','inferior'};
splineKnots = 15;
splineOrder = 4;

% Loop across the cardinal meridians and combine the RGC data. Use this
% aggregate data to define the knots for the spline fit.
aggregatePosition=[];
aggregateDensity=[];

for mm=1:length(cardinalMeridianAngles)

    mName = cardinalMeridianNames{mm};

    % Grab the density for this meridian
    coneDensitySqMm.(mName) = rawCellDensity.(mName);

    % Grab a copy of the support
    supportMm.(mName) = rawCellDensity.support;

    % handle leading zeros and trailing nans in the density vector
    coneDensitySqMm.(mName) = handleZerosAndNans(coneDensitySqMm.(mName));

    % handle nans around the optic disc in the density vector
    isvalididx=find(~isnan(coneDensitySqMm.(mName)));
    supportMm.(mName) = supportMm.(mName)(isvalididx);
    coneDensitySqMm.(mName) = coneDensitySqMm.(mName)(isvalididx);

    % make the initial support position slightly non-zero
    zeroIdx = find(supportMm.(mName)==0);
    if ~isempty(zeroIdx)
        supportMm.(mName)(zeroIdx) = 1e-6;
    end

    % aggregate the values across meridians
    aggregatePosition = [aggregatePosition supportMm.(mName)];
    aggregateDensity = [aggregateDensity coneDensitySqMm.(mName)];
end

% perform a least-squares spline fit with the specified knots and order
BformSplineFit=spap2(splineKnots, splineOrder, aggregatePosition', aggregateDensity');
knots = BformSplineFit.knots;

% Loop across the cardinal meridians again, and now perform the spline fit
% with the specified knots
for mm=1:length(cardinalMeridianAngles)
    mName = cardinalMeridianNames{mm};
    % Perform the spline fit
    BformSplineFit = ...
        spap2(knots, splineOrder, supportMm.(mName)', coneDensitySqMm.(mName)');
    % convert from B-form to piecewise polynomial form
    ppFormSplineFits{mm} = fn2fm(BformSplineFit,'pp');

end

% We will now assemble a new coefficient matrix for the spline fit that is
% an interpolation between the coefficient matrices obtained for the
% cardinal meridians. This is done by going through each element of the
% coefficient matrix and obtaining the four values, one from each meridian.
% The value from polar angle 0 is replicated and assigned a polar angle of
% 360; this allows the interpolation to wrap-around. We then fit a 'pchip'
% (Piecewise Cubic Hermite Interpolating Polynomial) to the set of 5 values
% and interpolate to the called-for, intermediate polar angle.

% make an empty coefficient matrix
interpCoefs = zeros(size(ppFormSplineFits{1}.coefs));

% replicate the values for polar angle zero at polar angle 360 to allow a
% wrap-around fit
angleBase=[cardinalMeridianAngles 360];

% loop over each element of the coefficient matrix
for cc=1:numel(interpCoefs)
    % obtain the set of coefficients from across the meridians
    thisCoefByMeridian = cellfun(@(x) x.coefs(cc), ppFormSplineFits);
    % duplicate the first coefficient (assuming it is polar angle zero)
    thisCoefByMeridian = [thisCoefByMeridian thisCoefByMeridian(1)];
    % fit the pchip model
    coefByAngleFit = fit(angleBase',thisCoefByMeridian','pchipinterp');
    % Store this function that will return coefficients for an arbitrary
    % polar angle
    interpCoefFunctions{cc} = coefByAngleFit;
end

% Create an anonymous function that returns a spline object that provides
% an interpolated piecewise polynomial fit structure as a function of polar
% angle
myCoefs = @(theta) reshape(cellfun(@(x) x(theta),interpCoefFunctions),size(interpCoefs));
mySpline = @(theta) setfield(ppFormSplineFits{1},'coefs',myCoefs(theta));

% Create an anonymous function using the interpolated spline. This function
% will return rgc density as a function of eccentricity in degrees. The
% transpose operations are needed so that that the function returns a row
% vector of density in response to a row vector of eccentricity support.

% Threshold y axis values below one to zero, and set density values close
% to the fovea to zero, as there is some weirdness in the spline at very
% low levels.
threhsoldY = @(y) y.*(y>1);
fitConeDensityCurcioMm = @(thetaDeg,eccentricityMM) threhsoldY(fnval(mySpline(thetaDeg),eccentricityMM')');

end


%%%% LOCAL FUNCTIONS

function densityVector = handleZerosAndNans(densityVector)
zeroIdx = find(densityVector(1:10)==0);
if ~isempty(zeroIdx)
    replacementVals = densityVector(zeroIdx(end)+1) ./ ((max(zeroIdx)-zeroIdx+1).^10.*10);
    if length(replacementVals)==1
        replacementVals = replacementVals./1e12;
    end
    densityVector(zeroIdx)=replacementVals;
end
nanIdx = find(isnan(densityVector(end-10:end)));
if ~isempty(nanIdx)
    densityVector(nanIdx+end-11)=1./(10.^(nanIdx-min(nanIdx)));
end
end

