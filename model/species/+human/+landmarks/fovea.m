function fovea = fovea( eye )
% Returns the fovea landmark sub-field of an eye model structure
%
% Syntax
%  fovea = human.landmarks.fovea( eye )
%
% Description:
%   Calculates the position on the retinal surface of the fovea. This
%   routine requires that the vertex subfield of landmarks has been
%   defined.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   fovea                 - Structure with the subfields degField,
%                           geodetic, and coords
%
% Examples:
%{
    % Plot the retinal surface and the positions of the fovea
    eye = modelEyeParameters('sphericalAmetropia',-1.5,'calcLandmarkFovea',true);
    S = eye.retina.S;
    boundingBox = eye.retina.boundingBox;
    figure
    quadric.plotSurface(S,boundingBox,[0.9 0.9 0.9],0.8,'b','g');
    camlight
    lighting gouraud
    hold on

    % Add the retinal landmarks
    plot3(eye.landmarks.vertex.coords(1),eye.landmarks.vertex.coords(2),eye.landmarks.vertex.coords(3),'+m','MarkerSize',10);
    plot3(eye.landmarks.fovea.coords(1),eye.landmarks.fovea.coords(2),eye.landmarks.fovea.coords(3),'+r','MarkerSize',10);

    % Add the geodetic path
    [geoDistance,~,~,geodeticPathCoords] = quadric.panouGeodesicDistance(S,eye.landmarks.fovea.geodetic,eye.landmarks.opticDisc.geodetic);
    plot3(geodeticPathCoords(:,1),geodeticPathCoords(:,2),geodeticPathCoords(:,3),'-y','MarkerSize',10);
%}
%{
    % Calculate the projective distance between the optic disc and fovea for
    % a range of spherical refractive errors
    SRvals = -10:1:2;
    for ii = 1:length(SRvals)
        eye = modelEyeParameters('sphericalAmetropia',SRvals(ii));
        odf(ii) = sqrt(sum((eye.landmarks.fovea.coords(2:3) - eye.landmarks.opticDisc.coords(2:3)).^2));
    end
    figure
    plot(SRvals,odf,'-*r');
    hold on
    % Compare these model values to the measurements reported by Jonas 2015
    % PloS One
    axialLengthFromSR = @(SR) 23.58-(SR.*0.299);
    odf_jonas = @(SR) 0.04 + 0.2.*(axialLengthFromSR(SR));
    plot(SRvals,odf_jonas(SRvals),'*b')
%}


% Check for the vertex subfield
if ~isfield(eye.landmarks,'vertex')
    error('human_landmarks_fovea:noVertexDefined','The vertex landmark must be defined in the eye structure');
end

% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% The calculation of the position of the retinal landmarks is based upon
% empirical measurements of visual field position. These measurements are
% done in the un-accommodated state, in the visible spectrum, with the eye
% and targets in air. It is also assumed that the position of the retinal
% landmarks were fixed at the point of maturity of the visual system.
% Therefore, the meta values of the passed eye structure are changed here
% to reflect these circumstances for this calculation.
eye.meta.spectralDomain = 'vis';
eye.meta.ageYears = 18;
eye.meta.accommodationDiopeters = 0;
cameraMedium = 'air';

% Set some fmincon options we will be using below
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% Here I assign the location of the fovea on the retina. This is done by
% first assigning an angle alpha between the optical axis of the eye and a
% second axis that originates in the fovea and arrives at a fixation point.
% This axis is practically realized here as the ray that connects the fovea
% through the center of the aperture stop and then exits the cornea at the
% angles alpha. This is in spirit similar to the "line of sight" axis,
% which has the property of passing through the center of the entrance
% pupil. The angle alpha is technically defined as the angle between the
% optical axis and the visual axis, which is the ray that originates from
% the fovea, passes through the nodal points of the eye, and arrives at the
% fixation point. As a practical matter, all three of these foveally-based
% axes are very similar to one another.
%
% The angle alpha has been proposed to vary with spherical refractive
% error:
%
%   Tabernero, Juan, et al. "Mechanism of compensation of aberrations in
%   the human eye." JOSA A 24.10 (2007): 3274-3283.
%
%   Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil shape as
%   viewed along the horizontal visual field." Journal of vision 13.6
%   (2013): 3-3.
%
% In each case, there is monotonic decrease in alpha with more negative
% spherical ametropia. The precise form of this is difficult to determine
% across studies. As a practical matter, I set the horizontal alpha value
% following the approach given in Figure 14 of Tabernero:
%
%{
    alpha0 = 5.8;
    L = @(SR) 16.5 / (16.5 - 0.299*SR );
    taberneroAlpha = @(SR) atand(L(SR)*tand(alpha0));
%}
%
% where SR is spherical refractive error in diopters and alpha0 is the
% alpha value for an emmetropic eye. The alpha0 value is 5.8, which is both
% the value given by Mathur 2013 and the median kappa value found in
% emmetropes in Hashemi 2010.
%{
    hashemiData = [6.75 5.75 4.38 4 3.88 3.75 3.38 3.13 3.13 2.88 2.63 2.63 2.5 2.38 2.25 2.25 2.25 2.13 2.13 1.88 1.88 1.75 1.75 1.75 1.63 1.63 1.63 1.63 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.38 1.38 1.38 1.38 1.38 1.38 1.38 1.38 1.25 1.25 1.25 1.25 1.25 1.25 1.13 1.13 1.13 1.13 1.13 1.13 1.13 1.13 1.13 1.13 1.13 1 1 1 1 1 1 1 1 1 1 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.88 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.38 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0 0 0 0 0 0 0 0 -0.13 -0.13 -0.13 -0.13 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.37 -0.38 -0.38 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.63 -0.63 -0.63 -0.75 -0.75 -0.75 -0.75 -0.75 -0.75 -0.75 -0.88 -0.88 -0.88 -0.88 -0.88 -0.88 -0.88 -0.88 -0.94 -1 -1 -1.13 -1.13 -1.25 -1.25 -1.25 -1.25 -1.38 -1.38 -1.5 -1.5 -1.5 -1.75 -1.88 -1.88 -2 -2.13 -2.63 -2.75 -2.75 -2.88 -3.13 -3.25 -3.25 -3.38 -3.75 -3.88 -3.88 -4.13 -4.13 -4.88 -5.25 -5.25 -5.38 -5.5 -6.13 -6.63 -8.13 -8.63 -9.25 -15.38; 8.29 3.09 5.38 6.83 7.15 4.69 4.06 4.1 7.76 4.87 5.05 3.36 6.12 5.11 7.14 4.58 5.33 5.75 5.3 6.85 4.79 4.08 3.66 6.55 5.5 4.82 5.02 4.83 5.09 7.66 3.84 6.09 6.4 4.99 6.31 3.22 5.81 7.11 4.24 5.54 2.1 6.88 4.52 4.18 6.96 6.58 5.56 6.34 6.18 4.86 5.43 5.82 6.14 6.48 6.06 6.11 4.17 5.29 6.68 6.08 6.28 7.13 5.85 7.04 3.7 7.31 5.54 5.11 4.78 6.04 5.9 8.06 6.88 3.32 3.32 4.14 5.54 4.68 5.13 5.36 6.41 5.29 5.63 3.04 3.96 5.92 5.73 5.3 6.19 4.42 4.94 6.49 5.72 4.74 4.57 3.46 6.07 4.91 5.2 5.5 4.5 6.98 5.73 10.75 6.79 5.82 2.71 3.49 7.84 5.9 4.78 2.71 4.93 6.52 5.95 5.21 5.23 4.87 5.59 6.04 6.77 5.67 5.47 5.26 4.49 3.78 6.7 5.36 6.14 7.4 5.14 5.79 4.04 5.95 6.69 4.7 7.52 6.73 5.83 5.81 5.43 6.65 6.49 5.63 5.24 3.97 5.96 6.19 7.67 6.26 7.69 5.83 6.05 5.47 6.41 3.72 5.41 4.43 7.02 5.51 4.31 5.69 6 6.07 4.76 8.27 5.23 4.56 5.85 7.39 6.15 4.9 6.1 5.63 6.13 6.49 6.19 6.31 5.34 6.36 4.42 6.86 5.98 6.63 3.19 4.56 4.81 9.1 4.96 5.57 4 5.8 6.26 4.87 6.49 3.58 4.38 5.71 6.86 4.98 5.18 5.51 6.86 5.5 5.79 5.97 7.06 5.45 6.5 4.96 6.27 4.43 4.78 7.09 5.66 4.88 3.5 6.18 3.49 6.24 5.6 5.69 5.69 7.18 4.81 5.87 5.88 7.65 5.27 7.15 4.66 5.95 5.88 5.74 6.99 5.49 4.01 7.02 4.55 6.07 4.2 4.88 4.54 6.62 7.71 5.57 5.48 1.11 5.11 6.43 6.68 6.19 3.73 7.88 4.49 5.42 5.18 5.92 4.75 6.03 4.58 6.61 6.04 2.94 2.35 4.06 3.48 3.86 4.16 5.4 6.65 3.92 5.2 3.3 3.97 0.94 2.1 5.08 2.37];
    taberneroData = [ 22.4046 22.6534 22.2886 22.2222 23.7313 24.6766 25.2736 24.4444 22.4378 22.2554 23.0680 23.0348 23.7479 24.5937 24.7761 25.4726 25.7877 27.8773; 8.8219 8.4110 7.5342 6.7123 7.1233 6.5205 6.1096 6.0274 5.6438 4.7945 4.1644 4.0274 3.9452 3.7260 4.8219 4.2192 3.7260 2.0274];
    % Convert axial length in taberneroData to SR
    ametropiaFromLength = @(x) (23.58 - x)./0.299;
    taberneroData(1,:) = ametropiaFromLength(taberneroData(1,:));
    % Obtain the mean emmetropic alpha
    fo = fitoptions('Method','NonlinearLeastSquares');
    k0=median(hashemiData(2,logical((hashemiData(1,:)>-0.5).*(hashemiData(1,:)<0.5))))
    expFunc = fittype( @(v,n,x) k0 .* ((x-v)./(-v)).^n,'independent','x','dependent','y','options',fo);
    dataX = [hashemiData(1,:) taberneroData(1,:)]';
    dataY = [hashemiData(2,:) taberneroData(2,:)]';
    expFit = fit(dataX,dataY,expFunc, ...
        'StartPoint',[16,0.5], ...
        'Lower',[-30,.5],'Upper',[-16,3],'Robust','Bisquare')
    plot(hashemiData(1,:),hashemiData(2,:),'or');
    hold on
    plot(taberneroData(1,:),taberneroData(2,:),'ob');
    plot(expFit.v:0.1:10,expFit(expFit.v:0.1:10),'-k')
%}
% For the vertical alpha, I assume an elevation of 3 degrees in the
% emmetropic eye.
%
a0 = [5.8 3.0 0];
L = @(SR) 16.5 / (16.5 - 0.299*SR );
alpha = @(SR) atand(L(SR).*tand(a0));
fovea.degField = alpha(eye.meta.sphericalAmetropia);
switch eye.meta.eyeLaterality
    case 'Right'
        % No change needed
    case 'Left'
        fovea.degField(1) = -fovea.degField(1);
    otherwise
        error('eye laterality not defined')
end

% Now calculate the location on the retina corresponding to this visual
% angle. The objective function is the difference in the visual field
% position of a candidate retinal point (G) and the desired value. Note
% that the elevational angle is inverted. This is because a negative value
% in this context corresponds to deflection of the visual axis upwards in
% the visual field.

% Some machinery to extract the second returned value from calcVisualAngle
theseAngles = @(G) wrapCalcVisualAngle(eye,eye.landmarks.vertex.geodetic,G,cameraMedium);
myObj = @(G) sum(angdiff(deg2rad(theseAngles(G)),deg2rad(fovea.degField(1:2).*[1 -1])).^2).*1e100;

% Define an x0 based upon laterality and quadric dimensions
switch eye.meta.eyeLaterality
    case 'Right'
        if isequal(quadric.dimensionSizeRank(S),[1 2 3]) || ...
                isequal(quadric.dimensionSizeRank(S),[1 3 2])
            x0 = [-85 -110 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[2 1 3])
            % Need to handle this case for extreme myopia
            error('The landmarks routine cannot model this extreme degree of myopia; exiting');
        end
    case 'Left'
        if isequal(quadric.dimensionSizeRank(S),[1 2 3]) || ...
                isequal(quadric.dimensionSizeRank(S),[1 3 2])
            x0 = [-85 110 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[2 1 3])
            % Need to handle this case for extreme myopia
            error('The landmarks routine cannot model this extreme degree of myopia; exiting');
        end
end

% Set the bounds
lb = [-89 -180 0];
ub = [-80 180 0];

% Perform the search
fovea.geodetic = fmincon(myObj, x0, [], [], [], [], lb, ub, [], opts);

% Obtain the visual axis ray
[~, ~, ~, outputRay1, ~, rayPath1] = calcVisualAngle(eye,eye.landmarks.vertex.geodetic,fovea.geodetic,[],[],cameraMedium);
fovea.outputRay = outputRay1;
fovea.rayPath = rayPath1;

% Obtain the coords coordinates of the fovea
fovea.coords = quadric.ellipsoidalGeoToCart(fovea.geodetic,S)';

end

%% LOCAL FUNCTION

% A miserable little function to return the second output of the
% calcVisualAngle function. Why oh why is this not built into MATLAB?
function visualAngleByPlane = wrapCalcVisualAngle(eye,G0,G1,cameraMedium)
[~, visualAngleByPlane ] = calcVisualAngle(eye,G0,G1,[],[],cameraMedium);
end