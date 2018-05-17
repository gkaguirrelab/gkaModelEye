function [transparentEllipseParams, RMSE, constraintError] = constrainedEllipseFit(Xp, Yp, lb, ub, nonlinconst)
% Non-linear fitting of an ellipse to a set of points
%
% Syntax:
%  [transparentEllipseParams, RMSE, constraintError] = constrainedEllipseFit(Xp, Yp, lb, ub, nonlinconst)
%
% Description:
%   The routine fits an ellipse to data by minimizing point-to-curve
%   distance, using an iterative procedure. The search is conducted over
%   the "transparent" ellipse parameterization, which has explicit values
%   for area eccentricity (aspect ratio), and theta. This allows us to set
%   boundaries and non-linear constraints upon these aspects of the fit.
%
%   This routine is a heavily modified version of a non-linear ellipse
%   fitting routine found within the "quadfit" matlab central toolbox. This
%   routine is dependent upon the quadfit toolbox.
%
% Input:
%	Xp, Yp                - nx1 Column vectors of the n points to be fit
%   lb, ub                - 1x5 vectors of the upper and lower bounds for
%                           the fit search, in transparent ellipse form
%   nonlinconst           - Function handle to a non-linear constraint
%                           function. This function should take as input
%                           the set of ellipse parameters in transparent
%                           form and return [c, ceq], where the optimizer
%                           constrains the solution such that c<=0 and
%                           ceq=0. This is an optional input; can be set
%                           as empty.
%
% Output:
%   transparentEllipseParams - Parameters of the best fitting ellipse
%                           expressed in transparent form [1x5 vector]
%   RMSE                  - Root mean squared error of the distance of each
%                           point in the data to the fitted ellipse
%   constraintError       - The value of the nonlinear constraint function
%                           for the best fitting ellipse
%
% Examples:
%{
    %% Calculate the time to perform a constrained fit
    % Generate some random transparent ellipse parameters
    nEllipses = 500;
    params=[(rand(nEllipses,1)-0.5)*1000, (rand(nEllipses,1)-0.5)*1000, 600+(rand(nEllipses,1)-0.5)*1000, rand(nEllipses,1)*0.75, rand(nEllipses,1)*pi];
    for ii = 1:nEllipses
        [xPoints(ii,:), yPoints(ii,:)] = ellipsePerimeterPoints(params(ii,:) );
    end
    tic
    for ii = 1:nEllipses
        fitParams(ii,:) = constrainedEllipseFit(xPoints(ii,:)', yPoints(ii,:)', [-500 -500 100 0 0], [500 500 1100 1 pi], [] );
    end
    timePerFitMsecs = toc / nEllipses * 1000;
    fprintf('Ellipse fitting time is %4.2f msecs.\n',timePerFitMsecs);
    e = max(abs(params-fitParams));
    fprintf('The maximum absolute fitting errors were [%4.3f %4.3f %4.3f %4.3f %4.3f].\n',e(1),e(2),e(3),e(4),e(5));
    fprintf('Note that elevated errors in theta fitting result from theta being unconstrained at low eccentricities.\n');
%}


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('ub',@isnumeric);
p.addRequired('lb',@isnumeric);
p.addRequired('nonlinconst',@(x) (isempty(x) || isa(x, 'function_handle')) );

% Parse and check the parameters
p.parse(Xp, Yp, ub, lb, nonlinconst);


%% Make an initial guess at the ellipse parameters
% This attempt is placed in a try-catch block, as the attempt can fail and
% return non-real numbers.
try
    % We sometimes obtain a singular matrix warning here; turn it
    % off temporarily
    warningState = warning;
    warning('off','MATLAB:singularMatrix');
    % use direct least squares ellipse fit to obtain an initial
    % estimate
    pInitImplicit = ellipsefit_direct(Xp,Yp);
    % Restore the warning state
    warning(warningState);
    % convert the initial estimate from implicit form to transparent form
    pInitTransparent = ellipse_ex2transparent(ellipse_im2ex(pInitImplicit));
    % place theta within the range of 0 to pi
    if pInitTransparent(5) < 0
        pInitTransparent(5) = pInitTransparent(5)+pi;
    end
catch
    % We couldn't find anything vaguely elliptical; return nans
    transparentEllipseParams=nan(1,5);
    RMSE=nan;
    constraintError=nan;
    return
end


%% Define the objective function
% This is the RMSE of the distance values of the boundary points to the
% ellipse fit
myFun = @(p) sqrt(nanmean(ellipsefit_distance(Xp,Yp,ellipse_transparent2ex(p)).^2));

% If the bounds and the nonlinear constraint function are all empty, then
% just return the initial estimate (and RMSE) obtained by direct fitting
if isempty(ub) && isempty(lb) && isempty(nonlinconst)
    transparentEllipseParams = pInitTransparent;
    RMSE = myFun(transparentEllipseParams);
    constraintError = nan;
    return
end


%% Perform non-linear search for transparent ellipse params

% Force the X and Y values of the initial guess to satisfy the nonlinear
% constraint
if ~isempty(nonlinconst)
    [~, ~, projectedEllipseOnImagePlane] = nonlinconst(pInitTransparent);
    pInitTransparent(4:5)=projectedEllipseOnImagePlane(4:5);
end

% define some search options
options = optimoptions(@fmincon,...
    'Algorithm','interior-point',...
    'Diagnostics','off',...
    'Display','off',...
    'DiffMinChange', 0.001);

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');

% Perform the non-linear search
[transparentEllipseParams, RMSE, ~, output] = ...
    fmincon(myFun, pInitTransparent, [], [], [], [], lb, ub, nonlinconst, options);

% Extract the constraint error from the output structure
constraintError = output.constrviolation;

% If we are close to zero for eccentricity, we may be in a local minimum
% that tends to find circles. Try a multiStart search with an increased
% lower bound on eccentricity.
if transparentEllipseParams(4) < 1e-12
    adjustedLB = lb;
    adjustedLB(4) = min([max([lb(4) 0.1]) ub(4)]);
    problem = createOptimProblem('fmincon','x0',pInitTransparent,...
        'objective',myFun,'lb',adjustedLB,'ub',ub,...
        'nonlcon',nonlinconst,'options',options);
    ms = MultiStart('Display','off','MaxTime',30,'StartPointsToRun','bounds-ineqs');
    [mstransparentEllipseParams, msRMSE, ~, msOutput] = run(ms,problem,20);
    if msRMSE <= RMSE
        RMSE = msRMSE;
        transparentEllipseParams = mstransparentEllipseParams;
        constraintError = msOutput.constrviolation;
    end
end

% Restore the warning state
warning(warningState);

end % function -- constrainedEllipseFit

