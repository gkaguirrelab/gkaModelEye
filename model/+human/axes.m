function axes = axes( eye, visualAxisDegRetina, opticDiscAxisDegRetina )

%% Axes - optical
% Eye axes are specified as rotations (in degrees) within the eye
% world coordinate frame for azimuth, elevation, and rotation. Axes
% are defined relative to the optical axis, which itself is set to
% be aligned with the p1 dimension of the eye world coordinate
% frame.
axes.optical.degRetina = [0 0 0];
axes.optical.mmRetina = [0 0 0];
axes.optical.degField = [0 0 0];


%% Axes - visual and blind spot
% The model establishes the position of the fovea and then sets the
% optic disc at a constant distance from the fovea in units of
% retinal degrees. The lines that connect these points on the fovea
% to the posterior nodal point of the eye define the visual and
% blind spot axes, respectively. The difference between these gives
% the position of the blind spot relative to fixation.
%
% Find the azimuthal arc in deg retina that produces a blind spot
% position in the horizontal and vertical directions that is equal
% to specified values from the literature. Values taken from Safren
% 1993 for their dim stimulus, under the assumption that this will
% be the most accurate given the minimization of light scatter. We
% model the fovea as being 3x closer to the optical axis than is
% the optic disc.
%{
    % Position of the blind spot in degrees of visual field
    % relative to fixation
    targetBlindSpotAngle = [-16.02 -1.84 0];
    blindSpotAngle = @(eye) axes.opticDisc.degField - axes.visual.degField;
    myObj = @(x) sum((blindSpotAngle(modelEyeParameters('opticDiscAxisDegRetina',[3/4*x(1),x(2)/2,0],'visualAxisDegRetina',-[1/4*x(1),x(2)/2,0])) - targetBlindSpotAngle).^2);
    options = optimoptions('fmincon','Display','off');
    retinalArcDeg = fmincon(myObj,[20 4],[],[],[],[],[],[],[],options);
    fprintf('Distance between the fovea and the center of the optic disc in retinal degrees in the right eye:\n');
    fprintf('\tazimuth = %4.4f; elevation = %4.4f \n\n', retinalArcDeg([1 2]));
%}
switch eye.meta.eyeLaterality
    case 'Right'
        opticDisc_WRT_foveaDegRetina = [-22.9384, 2.6078 ,0];
    case 'Left'
        opticDisc_WRT_foveaDegRetina = [22.9384, 2.6078 ,0];
end

% We next require the position of the fovea with respect to the
% optic axis in the emmetropic eye. We identify the position (in
% retinal degrees) of the fovea that results in a visual axis that
% has resulting alpha angles that match empirical results.  We
% assume an azimuth alpha of 5.8 degrees for an emmetropic eye
% (Figure 8 of Mathur 2013). We assume an elevation alpha of 2.5
% degrees, as this value, when adjusted to account for the longer
% axial length of the subjects in the Mathur study, best fits the
% Mathur data. Given these angles, we then calculate the
% corresponding position of the fovea w.r.t. the the optical axis
% of the eye (adjusted for eye laterality).
%{
    eye = modelEyeParameters();
    % These are the visual axis angles for an emmetropic eye
    targetAlphaAngle = [5.8  2.5  0];
    myComputedAlphaAzi = @(eye) axes.visual.degField(1);
    myObj = @(x) (targetAlphaAngle(1) - myComputedAlphaAzi(modelEyeParameters('visualAxisDegRetina',[x 0 0])))^2;
    aziFoveaEmmetropic = fminsearch(myObj,9)
    myComputedAlphaEle = @(eye) axes.visual.degField(2);
    myObj = @(x) (targetAlphaAngle(2) - myComputedAlphaEle(modelEyeParameters('visualAxisDegRetina',[aziFoveaEmmetropic x 0])))^2;
    eleFoveaEmmetropic = fminsearch(myObj,2)
%}
switch eye.meta.eyeLaterality
    case 'Right'
        fovea_WRT_opticAxisDegRetina_emmetrope = [8.2964 -3.5762 0];
    case 'Left'
        fovea_WRT_opticAxisDegRetina_emmetrope = [-8.2964 -3.5762 0];
end

% In our model, the fovea moves towards the apex of the posterior
% chamber as the eye becomes closer to spherical. We implement this
% effect by calculating the ratio of the vitreous chamber axes.
%{
    format long
    probeEye = modelEyeParameters('sphericalAmetropia',0);
    eccen_p1p2 = (1-probeEye.vitreousChamber.radii(1)/probeEye.vitreousChamber.radii(2))
    eccen_p1p3 = (1-probeEye.vitreousChamber.radii(1)/probeEye.vitreousChamber.radii(3))
    format
%}
foveaPostionScaler(1) = (1-eye.vitreousChamber.radii(1)/eye.vitreousChamber.radii(2))/0.111716335829885;
foveaPostionScaler(2) = (1-eye.vitreousChamber.radii(1)/eye.vitreousChamber.radii(3))/0.105571718627770;
foveaPostionScaler(3) = 1;
axes.visual.degRetina = fovea_WRT_opticAxisDegRetina_emmetrope.*foveaPostionScaler;

% The optic disc maintains a fixed distance (in retinal degrees)
% from the fovea
axes.opticDisc.degRetina = opticDisc_WRT_foveaDegRetina + axes.visual.degRetina;

% If a visualAxisDegRetina or opticDiscAxisDegRetina key-value pair
% was passed, override the computed value. This is used primarily
% during model development.
if ~isempty(visualAxisDegRetina)
    axes.visual.degRetina = visualAxisDegRetina;
end
if ~isempty(opticDiscAxisDegRetina)
    axes.opticDisc.degRetina = opticDiscAxisDegRetina;
end

% Calculate the foveal and optic disc positions in terms of mm of
% retina. This requires the elliptic integral. The parameter
% "theta" has a value of zero at the apex of the ellipse along the
% axial dimension (p1).
ellipticIntegral_p1p2=@(theta) sqrt(1-sqrt(1-eye.vitreousChamber.radii(2).^2/eye.vitreousChamber.radii(1).^2)^2.*(sin(theta)).^2);
ellipticIntegral_p1p3=@(theta) sqrt(1-sqrt(1-eye.vitreousChamber.radii(3).^2/eye.vitreousChamber.radii(1).^2)^2.*(sin(theta)).^2);
arcLength_p1p2 = @(theta1,theta2) eye.vitreousChamber.radii(1).*integral(ellipticIntegral_p1p2, theta1, theta2);
arcLength_p1p3 = @(theta1,theta2) eye.vitreousChamber.radii(1).*integral(ellipticIntegral_p1p3, theta1, theta2);

% For the calculation, the first theta value is zero, as we are
% calculating distance from the vitreous chamber apex (i.e., the
% intersection of the optical axis with the retina).
axes.visual.mmRetina = [arcLength_p1p2(0,deg2rad(axes.visual.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.visual.degRetina(2))), 0];
axes.opticDisc.mmRetina = [arcLength_p1p2(0,deg2rad(axes.opticDisc.degRetina(1))), arcLength_p1p3(0,deg2rad(axes.opticDisc.degRetina(2))), 0];

% Calculate the foveal position in eyeWorld coordinates.
phi = -axes.visual.degRetina(1);
theta = -axes.visual.degRetina(2);
x = eye.vitreousChamber.radii(1) * cosd(theta) * cosd(phi);
y = eye.vitreousChamber.radii(2) * cosd(theta) * sind(phi);
z = eye.vitreousChamber.radii(3) * sind(theta);

% Note this location in the vitreous chamber field
axes.visual.coords = [-x y -z] + eye.vitreousChamber.center;

% Calculate the optic disc position in eyeWorld coordinates.
phi = -axes.opticDisc.degRetina(1);
theta = -axes.opticDisc.degRetina(2);
x = eye.vitreousChamber.radii(1) * cosd(theta) * cosd(phi);
y = eye.vitreousChamber.radii(2) * cosd(theta) * sind(phi);
z = eye.vitreousChamber.radii(3) * sind(theta);

% Store this location
axes.opticDisc.coords = [-x y -z] + eye.vitreousChamber.center;

% Calcuate the optic disc and visual axes in deg of visual field,
% using the nodal point of the eye. For the visual axis, these
% values correspond to alpha / kappa, the angles between the visual
% and optical /pupillary axes. The difference between the visual
% and optic disc axes specifies the location of the physiologic
% blind spot relative to fixation.
axes.visual.degField(1) = atand((axes.visual.coords(2) - eye.lens.nodalPoint(2)) / (axes.visual.coords(1) - eye.lens.nodalPoint(1)));
axes.visual.degField(2) = -(-atand((axes.visual.coords(3) - eye.lens.nodalPoint(3)) / (axes.visual.coords(1) - eye.lens.nodalPoint(1))));
axes.visual.degField(3) = 0;
axes.opticDisc.degField(1) = atand((axes.opticDisc.coords(2) - eye.lens.nodalPoint(2)) / (axes.opticDisc.coords(1) - eye.lens.nodalPoint(1)));
axes.opticDisc.degField(2) = -(-atand((axes.opticDisc.coords(3) - eye.lens.nodalPoint(3)) / (axes.opticDisc.coords(1) - eye.lens.nodalPoint(1))));
axes.opticDisc.degField(3) = 0;

end

