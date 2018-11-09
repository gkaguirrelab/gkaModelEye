function [visualAngles, rayPath0, rayPath1, totalAngle ] = calcVisualAngle(eye,G0,G1,X0,X1,cameraMedium)
% The visual angles between two retinal points
%
% Syntax:
%  visualAngles = calcVisualAngle(sceneGeometry,G0,G1,X0,X1)
%
% Description
%   Given a sceneGeometry and two coordinates on the retinal surface, the
%   routine returns a vector that contains the visual angle (in degrees)
%   between the two points, projected on the p1p2 and p1p3 planes (i.e.,
%   horizontal and vertical visual angle).
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
%   visualAngles          - 1x2 vector with the visual angles, in degrees
%                           between the two points within the p1p2 and p1p3
%                           planes.
%   rayPath0, rayPath1    - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%
% Examples:
%{
    % Display a map of visual angle on the retinal surface
    eye = modelEyeParameters('eyeLaterality','os','skipNodalPoint',true);
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
    for beta = -90:3:0
        for omega = -180:5:180
            visualAngles = calcVisualAngle(eye,eye.axes.visual.geodetic,[beta omega 0]);
            eccen = sqrt(sum(visualAngles.^2));
            if ~isnan(eccen)
                if eccen > 90
                    colorTriple = [1 1 1];
                    markerSize = 25;
                else
                    colorTriple = c(round((eccen./90)*(nColors-1)+1),:);
                    markerSize = 50;
                end
                coord = quadric.ellipsoidalGeoToCart( [beta omega 0], S );
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
    for SR = -5:1:2
        eye = modelEyeParameters('sphericalAmetropia',SR,'skipNodalPoint',true);
        S = eye.retina.S;
        G0 = eye.axes.visual.geodetic;
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
[R0, rayPath0] = calcNodalRay(eye,[],X0,cameraMedium);
[R1, rayPath1] = calcNodalRay(eye,[],X1,cameraMedium);

% Calculate and return the signed angles between the two rays
[totalAngle, angle_p1p2, angle_p1p3] = quadric.angleRays( R0, R1 );
visualAngles = [angle_p1p2, angle_p1p3];

end