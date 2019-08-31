function [G,X,angleError] = calcRetinaFieldPoint( eye, degField, cameraMedium )
% Retinal point at a specified visual field location w.r.t. optical axis
%
% Syntax
%  [G,X,angleError] = calcRetinaFieldPoint( eye, degField, cameraMedium )
%
% Description:
%   Calculates the position on the retina of a point that has the specified
%   visual field position, where visual field position is defined with
%   respect to the optical axis of the eye. The identified retinal point
%   has the property that a nodal ray that emerges from this point exits
%   the cornea at angles with respect to the optical axis that are equal to
%   the values specified in degField.
%
% Inputs
%   eye                   - Structure.
%   degField              - 2x1 vector that specifies the visual field
%                           position of a point in degree angles along the
%                           horizontal and vertical meridian with respect
%                           to the optical axis of the eye.
%   cameraMedium          - Char vector. Used to determine the refractive
%                           index of the medium of the visual field.
%                           Defaults to "air" if not passed.
%
% Outputs
%   G                     - 3x1 vector that provides the geodetic
%                           coordinates beta, omega, and elevation in units
%                           of degrees. Beta is defined over the range
%                           -90:90, and omega over the range -180:180.
%                           Elevation has an obligatory value of zero as
%                           this solution is only defined on the surface.
%   X                     - 3x1 vector that specifies the Cartesian
%                           location of a point on the quadric surface.
%   angleError            - Scalar. The angle between the desired output
%                           ray and the solution. 
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    degField = [5.8 3.0];
    [G,X,angleError] = calcRetinaFieldPoint( sceneGeometry.eye, degField);
    [outputRay,rayPath] = calcNodalRay(sceneGeometry.eye,G);
    plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
%}
%{
    % Relation between axial length and mm per deg at the retinal apex
    deltaAngles=[sqrt(1/2)/2 sqrt(1/2)/2 0];
    length = [];
    mmPerDeg = [];
    for SR = -7:1:3
        eye = modelEyeParameters('sphericalAmetropia',SR);
        [~,X0] = calcRetinaFieldPoint( eye, -deltaAngles);
        [~,X1] = calcRetinaFieldPoint( eye, deltaAngles);
        length(end+1) = eye.meta.axialLength;
        mmPerDeg(end+1) = sqrt(sum((X0-X1).^2)) ./ sqrt(sum((deltaAngles.*2).^2));
    end
    figure
    plot(length,mmPerDeg,'xr');
    xlabel('axial length [mm]');
    ylabel('mm retina per deg visual angle');
    vals = polyfit(length,mmPerDeg,1);
    fprintf('Retinal mm per deg visual field at the viterous chamber apex = (%2.3f * axialLength) %2.3f \n',vals(1),vals(2));
%}

% Handle incomplete inputs
if nargin<2
    error('Need to specify an eye structure and the visual field angles');
end

if nargin==2
    cameraMedium = 'air';
end

% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% To set an x0 guess, identify the retinal point that is at the specified
% angle w.r.t. to the optical axis and the center of the aperture stop.
R = quadric.anglesToRay(eye.stop.center',degField(1),degField(2));
x0 = quadric.intersectRay(S,R,eye.retina.side,eye.retina.boundingBox);

% Convert the x0 zero guess into ellipsoidal geodetic coordinates, beta
% (latitude), omega (longitude), and elevation.
g0 = quadric.cartToEllipsoidalGeo( x0, S );

% Define the optical axis ray
opticalAxis = [0 1; 0 0; 0 0];

% The objective is to match the angles between the output nodal ray and the
% optical axis to the specified degField angles.  Note that the elevational
% angle is inverted. This is because a negative value in this context
% corresponds to deflection of the visual axis upwards in the visual field.
% We force the third component of the geodetic coordinate (elevation) to be
% zero.
myObj = @(G) sqrt(sum((-degField(1:2).*[1 -1]-wrapAngleRays(calcNodalRay(eye,[G; 0],[],cameraMedium),opticalAxis)).^2));

% define some search options
options = optimset;

% Perform the search
G = [0; 0; 0];
[G(1:2),angleError] = fminsearch(myObj, g0(1:2), options);

% Obtain the Cartesian coordinates of the fovea
X = quadric.ellipsoidalGeoToCart(G,S);

end

%% LOCAL FUNCTION

% Wrap 
function angles = wrapAngleRays(R0,R1)
[~, angles(1), angles(2) ] = quadric.angleRays(R0,R1);
end
