function visualAngles = calcVisualAngle(eye,G0,G1,X0,X1,cameraMedium)
% The visual angles between two retinal points
%
% Syntax:
%  visualAngle = calcVisualAngle(sceneGeometry,G0,G1,X0,X1)
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
%   visualAngles          - 1x2 vector with the visual angle, in degrees
%                           between the two points within the p1p2 and p1p3
%                           planes.
%
% Examples:
%{
    % Display a map of visual angle on the retinal surface
    eye = modelEyeParameters('eyeLaterality','os');
    S = eye.retina.S;
    boundingBox = eye.retina.boundingBox;
    figure
    quadric.plotSurface(S,boundingBox,[0.9 0.9 0.9],0.8,'b','g');
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
                else
                    colorTriple = c(round((eccen./90)*(nColors-1)+1),:);
                end
                plot3(coord(1),coord(2),coord(3),'.','MarkerSize',50,'Color',colorTriple);
            end
        end
    end
    % Add the retinal landmarks
    plot3(eye.axes.optical.coords(1),eye.axes.optical.coords(2),eye.axes.optical.coords(3),'+k','MarkerSize',10);
    plot3(eye.axes.visual.coords(1),eye.axes.visual.coords(2),eye.axes.visual.coords(3),'*k','MarkerSize',10);
    plot3(eye.axes.opticDisc.coords(1),eye.axes.opticDisc.coords(2),eye.axes.opticDisc.coords(3),'ok','MarkerSize',10);
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

% Handle to the virtual image function; use the MEX version if available
if exist('virtualImageFuncMex')==3
    refractionHandle = @virtualImageFuncMex;
else
    refractionHandle = @virtualImageFunc;
end

% Assemble arguments for the virtual image function to trace from the
% retina to the center of the pupil aperture. The "camera" position is set
% as the pupil center, with the dimensions rearranged for world
% coordinates.
args = {eye.pupil.center([2 3 1])', ...
    eye.rotationCenters, ...
    assembleOpticalSystem( eye, 'surfaceSetName','retinaToPupil', 'cameraMedium', cameraMedium )};

% Set eyePose to all zeros (i.e., the eye is not rotated).
eyePose = [0 0 0 0];

% Ray trace to the center of the pupil
R0 = refractionHandle(X0, eyePose, args{:});
R1 = refractionHandle(X1, eyePose, args{:});

% Normalize the rays, and reverse the direction so that the ray is headed
% from the pupil center out towards the cornea
R0 = quadric.normalizeRay(R0');
R0(:,2)=-R0(:,2);
R1 = quadric.normalizeRay(R1');
R1(:,2)=-R1(:,2);

% Ray trace from the pupil center through the cornea
R0 = rayTraceQuadrics(R0, assembleOpticalSystem( eye, 'surfaceSetName','pupilToCamera','cameraMedium',cameraMedium ));
R1 = rayTraceQuadrics(R1, assembleOpticalSystem( eye, 'surfaceSetName','pupilToCamera','cameraMedium',cameraMedium ));

% Calculate and return the signed angles between the two rays
[~, angle_p1p2, angle_p1p3] = quadric.angleRays( R0, R1 );
visualAngles = [angle_p1p2, angle_p1p3];

end