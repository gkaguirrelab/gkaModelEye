function [visualAngleTotal, visualAngleByPlane, outputRay0, outputRay1, rayPath0, rayPath1 ] = calcVisualAngle(eye,G0,G1,X0,X1,cameraMedium)
% The visual angles between two retinal points
%
% Syntax:
%  visualAngleTotal = calcVisualAngle(sceneGeometry,G0,G1,X0,X1)
%
% Description
%   Given a sceneGeometry and two coordinates on the retinal surface, the
%   routine returns the total visual angle (in degrees) between the two
%   points, as well as a vector that contains the angles projected on the
%   p1p2 and p1p3 planes (i.e., horizontal and vertical visual angle).
%
%   The routine can accept points on the ellipsoidal surface specified in
%   either Cartesian or ellipsoidal geodetic coordinates.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%   G0, G1                - 3x1 vectors that provide the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X0, X1                - 3x1 vectors that specify the Cartesian
%                           location of points on the quadric surface.
%   cameraMedium          - The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   visualAngleTotal      - Scalar. The angle in degrees of visual field
%                           between two points on the retina.
%   visualAngleByPlane    - 1x2 vector with the visual angles, in degrees
%                           between the two points within the p1p2 and p1p3
%                           planes.
%   outputRay0, outputRay1 - 3x2 matricies that specify the ray as a unit 
%                           vector of the form [p; d], corresponding to
%                               R = p + t*u
%                           where p is vector origin, d is the direction
%                           expressed as a unit step, and t is unity.
%   rayPath0, rayPath1    - 3xm matricies that provide the ray coordinates
%                           at each surface. The value for (e.g.)
%                           rayPath0(:,1) is equal to initial position. If
%                           a surface is missed, then the coordinates for
%                           that surface will be nan.
%
% Examples:
%{
    % Demonstrate the basic calculation using Cartesian coordiantes
    eye = modelEyeParameters('calcLandmarkFovea',true);
    visualAngleTotal = calcVisualAngle(eye, [], [], eye.landmarks.vertex.coords, eye.landmarks.fovea.coords);
%}
%{
    % Display a map of visual angle on the retinal surface
    eye = modelEyeParameters('eyeLaterality','os','calcLandmarkFovea',true);
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
    plotRezGeodeticDeg = [15, 20];
    for beta = -90:plotRezGeodeticDeg(1):0
        for omega = -180:plotRezGeodeticDeg(2):180
            visualAngleTotal = calcVisualAngle(eye,eye.landmarks.fovea.geodetic,[beta omega 0]);
            if ~isnan(visualAngleTotal)
                if visualAngleTotal > 90
                    colorTriple = [1 1 1];
                    markerSize = 25;
                else
                    colorTriple = c(round((visualAngleTotal./90)*(nColors-1)+1),:);
                    markerSize = 50;
                end
                coord = quadric.ellipsoidalGeoToCart( [beta omega 0], S );
                plot3(coord(1),coord(2),coord(3),'.','MarkerSize',markerSize,'Color',colorTriple);
            end
        end
    end
    % Add the retinal landmarks
    plot3(eye.landmarks.vertex.coords(1),eye.landmarks.vertex.coords(2),eye.landmarks.vertex.coords(3),'+k','MarkerSize',10);
    plot3(eye.landmarks.fovea.coords(1),eye.landmarks.fovea.coords(2),eye.landmarks.fovea.coords(3),'*k','MarkerSize',10);
%}
%{
    % Calculate deg/mm at the fovea for an emmetropic eye
    % length
    eye = modelEyeParameters('calcLandmarkFovea',true);
    S = eye.retina.S;
    G0 = eye.landmarks.fovea.geodetic;
    G1 = G0 + [0.1 0.1 0];
    visualAngleTotal = calcVisualAngle(eye,G0,G1);
    X0 = quadric.ellipsoidalGeoToCart(G0,S);
    X1 = quadric.ellipsoidalGeoToCart(G1,S);
    d = sqrt(sum((X0-X1).^2));
    mmPerDeg = d/visualAngleTotal;
    outline = sprintf('%2.2f mm/deg at the fovea for an emmetropic eye.\n',mmPerDeg);
    fprintf(outline);
%}
%{
    % Calculate deg/mm at the fovea as a function of ametropia and axial
    % length
    mmPerDeg = [];
    axialLengths = [];
    for SR = -5:1:2
        eye = modelEyeParameters('sphericalAmetropia',SR,'calcLandmarkFovea',true);
        S = eye.retina.S;
        G0 = eye.landmarks.fovea.geodetic;
        G1 = G0 + [0.1 0.1 0];
        visualAngleTotal = calcVisualAngle(eye,G0,G1);
        X0 = quadric.ellipsoidalGeoToCart(G0,S);
        X1 = quadric.ellipsoidalGeoToCart(G1,S);
        d = sqrt(sum((X0-X1).^2));
        mmPerDeg(end+1) = d/visualAngleTotal;
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


% Handle incomplete inputs
if nargin<=2
    error('Invalid number of input arguments');
end

if nargin==3
    % If only three input values were passed, derive the X0/X1 Cartesian
    % coordinates from the ellipsoidal geodetic coordinates.
    S = eye.retina.S;
    X0 = quadric.ellipsoidalGeoToCart( G0, S )';
    X1 = quadric.ellipsoidalGeoToCart( G1, S )';
    cameraMedium = 'air';
end

if nargin==4
    error('Invalid number of input arguments');
end

if nargin==5
    % Check if the X0/X1 values are empty. If so, derive the X0/X1
    % Cartesian coordinates from the ellipsoidal geodetic coordinates.
    if isempty(X0)
        S = eye.retina.S;
        X0 = quadric.ellipsoidalGeoToCart( G0, S )';
        X1 = quadric.ellipsoidalGeoToCart( G1, S )';
    end
    cameraMedium = 'air';
end

if nargin==6
    % Check if the X0/X1 values are empty. If so, derive the X0/X1
    % Cartesian coordinates from the ellipsoidal geodetic coordinates.
    if isempty(X0)
        S = eye.retina.S;
        X0 = quadric.ellipsoidalGeoToCart( G0, S )';
        X1 = quadric.ellipsoidalGeoToCart( G1, S )';
    end
end

% Obtain the output ray segments corresponding to the nodal ray for each
% retinal coordinate
[outputRay0, rayPath0] = calcNodalRay(eye,[],X0,cameraMedium);
[outputRay1, rayPath1] = calcNodalRay(eye,[],X1,cameraMedium);

% Calculate and return the signed angles between the two rays
[visualAngleTotal, angle_p1p2, angle_p1p3] = quadric.angleRays( outputRay0, outputRay1 );
visualAngleByPlane = [angle_p1p2, angle_p1p3];

end