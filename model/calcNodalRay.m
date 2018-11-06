function rayPath = calcNodalRay(eye,G0,X0,cameraMedium)
% Returns the path of the nodal ray from a retinal point
%
% Syntax:
%  rayPath = calcNodalRay(eye,G0,X0,cameraMedium)
%
% Description
%   Given a sceneGeometry and a coordinate on the retinal surface, the
%   routine returns a matrix that contains the path of a ray that satisfies
%   the property that the angle (wrt the optical axis) of the ray as it
%   departs the retina is equal to the angle with which it departs the
%   cornea. This is a "nodal ray":
%
%       Harris, W. F. "Nodes and nodal points and lines in eyes and other
%       optical systems." Ophthalmic and Physiological Optics 30.1 (2010):
%       24-42.
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   G0                    - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0                    - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   visualAngles          - 1x2 vector with the visual angles, in degrees,
%                           between the retinal point and the optical axis
%                           within the p1p2 and p1p3 planes.
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%
% Examples:
%{
    % Display a map of visual angle on the retinal surface
    eye = modelEyeParameters('eyeLaterality','os');
    S = eye.retina.S;
    boundingBox = eye.retina.boundingBox;
    figure
    quadric.plotSurface(S,boundingBox,[0.9 0.9 0.9],0.8);
    camlight
    lighting gouraud
    hold on
    % Sweep over ellipsoidal geodetic coordinates on the retinal surface
    % and obtain the visual angle and each point w.r.t. the fovea.
    c = jet();
    nColors = size(c,1);
    for beta = -90:3:40
        for omega = -180:5:180
            coord = quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
            visualAngles = calcVisualAngle(eye,eye.axes.visual.geodetic,[beta omega 0]);
            eccen = sqrt(sum(visualAngles.^2));
            if ~isnan(eccen)
                if eccen > 90
                    colorTriple = [1 1 1];
                    markerSize = 50;
                else
                    colorTriple = c(round((eccen./90)*(nColors-1)+1),:);
                    markerSize = 100;
                end
                plot3(coord(1),coord(2),coord(3),'.','MarkerSize',markerSize,'Color',colorTriple);
            end
        end
    end
    % Add the retinal landmarks
    plot3(eye.axes.optical.coords(1),eye.axes.optical.coords(2),eye.axes.optical.coords(3),'+k','MarkerSize',10);
    plot3(eye.axes.visual.coords(1),eye.axes.visual.coords(2),eye.axes.visual.coords(3),'*k','MarkerSize',10);
    plot3(eye.axes.opticDisc.coords(1),eye.axes.opticDisc.coords(2),eye.axes.opticDisc.coords(3),'ok','MarkerSize',10);
%}
%{
    % Calculate deg/mm at the fovea as a function of ametropia and axial
    % length
    mmPerDeg = [];
    axialLengths = [];
    for SR = -10:1:2
        eye = modelEyeParameters('sphericalAmetropia',SR);
        S = eye.retina.S;
        G0 = eye.axes.opticDisc.geodetic;
        G1 = G0 + [0.1 0.1 0];
        [~, ~, ~, totalAngle ] = calcVisualAngle(eye,G0,G1);
        X0 = quadric.ellipsoidalGeoToCart(G0,S);
        X1 = quadric.ellipsoidalGeoToCart(G1,S);
        d = sqrt(sum((X0-X1).^2));
        mmPerDeg(end+1) = d/totalAngle;
        radii = quadric.radii(eye.retina.S);
        center = quadric.center(eye.retina.S);
        axialLengths(end+1) = radii(1)-center(1);
    end
    figure
    plot(axialLengths,mmPerDeg,'or');
    p = polyfit(axialLengths,mmPerDeg,1);
    hold on
    x = linspace(22,28);
    y = polyval(p,x);
    plot(x,y,'-k');
    xlabel('axial length [mm]');
    ylabel('mm / degree visual angle at the fovea');
%}


if nargin<2
    error('Invalid number of input arguments');
end

if nargin==2
    % If only two input values were passed, derive the X0 Cartesian
    % coordinates from the ellipsoidal geodetic coordinates.
    S = eye.retina.S;
    X0 = quadric.ellipsoidalGeoToCart( G0, S )';
    cameraMedium = 'air';
end

if nargin==3
    % Check if the X0 value is empty. If so, derive the X0
    % Cartesian coordinates from the ellipsoidal geodetic coordinate.
    if isempty(X0)
        S = eye.retina.S;
        X0 = quadric.ellipsoidalGeoToCart( G0, S )';
    end
    cameraMedium = 'air';
end

if nargin==4
    % Check if the X0 value is empty. If so, derive the X0
    % Cartesian coordinates from the ellipsoidal geodetic coordinate.
    if isempty(X0)
        S = eye.retina.S;
        X0 = quadric.ellipsoidalGeoToCart( G0, S )';
    end
end

% Handle to the virtual image function; use the MEX version if available
if exist('virtualImageFuncMex')==3
    refractionHandle = @virtualImageFuncMex;
else
    refractionHandle = @virtualImageFunc;
end

% Assemble the optical system
opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName','retinaToCamera', 'cameraMedium', cameraMedium );

% Define some options for the fmincon call in the loop
options = optimoptions(@fmincon,...
    'Display','off');

        % Define an error function that reflects the difference in angles
        % from the initial ray and the output ray from the optical system
        myError = @(x) calcOffsetFromParallel(opticalSystem,assembleInputRay(coord,X0(1),X0(2)));

        % Supply an x0 guess as the ray that connects the retinal point
        % with the center of the pupil
        [~, angle_p1p2, angle_p1p3] = quadric.angleRays( [0 0 0; 1 0 0]', quadric.normalizeRay([coord'; eye.pupil.center-coord']') );
        angle_p1p2 = deg2rad(angle_p1p2);
        angle_p1p3 = -deg2rad(angle_p1p3);

        % Perform the search
        inputRayAngles = fmincon(myError,[angle_p1p2 angle_p1p3],[],[],[],[],[-pi/2,-pi/2],[pi/2,pi/2],[],options);

        % Calculate and save the outputRay and the raypath
        [outputRay,rayPath] = rayTraceQuadrics(assembleInputRay(coord,inputRayAngles(1),inputRayAngles(2)), opticalSystem);


% Calculate and return the signed angles between the two rays
[totalAngle, angle_p1p2, angle_p1p3] = quadric.angleRays( R0, R1 );
visualAngles = [angle_p1p2, angle_p1p3];

end


% Local function. Converts angles relative to the optical axis to a unit
% vector ray.
function inputRay = assembleInputRay(p,angle_p1p2,angle_p1p3)
u = [1; tan(angle_p1p2); tan(angle_p1p3)];
u = u./sqrt(sum(u.^2));
inputRay = [p, u];
end

% Local function. Performs the ray trace through the optical system of the
% eye and then calculates the angle between the initial ray and the output
% ray.
function angleError = calcOffsetFromParallel(opticalSystem,inputRay)
    exitRay = rayTraceQuadrics(inputRay, opticalSystem);
    angleError = quadric.angleRays( inputRay, exitRay );
end
