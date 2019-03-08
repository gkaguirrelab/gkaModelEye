function opticDisc = opticDisc( eye )
% Returns the opticDisc landmark sub-field of an eye model structure
%
% Syntax
%  opticDisc = human.landmarks.opticDisc( eye )
%
% Description:
%   Calculates the position on the retinal surface of the optic disc. This
%   routine requires that the vertex and fovea subfields of landmarks has
%   been defined.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   fovea                 - Structure with the subfields degField,
%                           geodetic, and coords
%



% Check for the vertex subfield
if ~isfield(eye.landmarks,'vertex')
    error('human_landmarks_opticDisc:noVertexDefined','The vertex landmark must be defined in the eye structure');
end

if ~isfield(eye.landmarks,'fovea')
    error('human_landmarks_opticDisc:noFoveaDefined','The fovea landmark must be defined in the eye structure');
end


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

% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% Set some fmincon options we will be using below
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

% The center of the optic disc is nasal and superior to the fovea.
% Rohrschneider measured the distance (in degrees of visual angle) between
% the center of the optic disc and the fovea in 104 people, and found a
% mean horizontal value of 15.5 +- 1.1 degrees SD, and a vertical value of
% 1.5 +- 0.9 degrees SD. There was no dependence of the value upon
% ametropia.
%
%   Rohrschneider, Klaus. "Determination of the location of the fovea on
%   the fundus." Investigative ophthalmology & visual science 45.9 (2004):
%   3257-3258.
%
% By adding the position of fovea in degrees relative to the optical axis,
% we obtain the position of the optic disc in degrees relative to the
% optical axis.
%
switch eye.meta.eyeLaterality
    case 'Right'
        opticDisc.degField = [-15.5 -1.5 0]+eye.landmarks.fovea.degField;
    case 'Left'
        opticDisc.degField = [15.5 -1.5 0]+eye.landmarks.fovea.degField;
end

% Now calculate the location on the retina corresponding to this visual
% angle. The objective function is the difference in the visual field
% position of a candidate retinal point (G) and the desired value. Note
% that the elevational angle is inverted. This is because a negative value
% in this context corresponds to deflection of the visual axis upwards in
% the visual field.

% Some machinery to extract the second returned value from calcVisualAngle
theseAngles = @(G) wrapCalcVisualAngle(eye,eye.landmarks.vertex.geodetic,G,cameraMedium);
myObj = @(G) sum(angdiff(deg2rad(theseAngles(G)),deg2rad(opticDisc.degField(1:2).*[1 -1])).^2).*1e100;

% Define an x0 based upon laterality and quadric dimensions
switch eye.meta.eyeLaterality
    case 'Right'
        if isequal(quadric.dimensionSizeRank(S),[1 2 3])
            x0 = [-72 90 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[1 3 2])
            x0 = [-88 -50 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[2 1 3])
            % Need to handle this case for extreme myopia
            error('Geoff needs to implement this case');
        end
    case 'Left'
        if isequal(quadric.dimensionSizeRank(S),[1 2 3])
            x0 = [-72 -90 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[1 3 2])
            x0 = [-88 50 0];
        end
        if isequal(quadric.dimensionSizeRank(S),[2 1 3])
            % Need to handle this case for extreme myopia
            error('Geoff needs to implement this case');
        end
end

% Set the search bounds
lb = [-89 -180 0];
ub = [-60 180 0];

% Perform the search
opticDisc.geodetic = fmincon(myObj, x0, [], [], [], [], lb, ub, [], opts);

% Obtain the coords of the optic disc
opticDisc.coords = quadric.ellipsoidalGeoToCart(opticDisc.geodetic,S)';


end



%% LOCAL FUNCTION

% A miserable little function to return the second output of the
% calcVisualAngle function. Why oh why is this not built into MATLAB?
function visualAngleByPlane = wrapCalcVisualAngle(eye,G0,G1,cameraMedium)
[~, visualAngleByPlane ] = calcVisualAngle(eye,G0,G1,[],[],cameraMedium);
end