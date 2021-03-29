function fovea = fovea( eye )
% Returns the fovea landmark sub-field of an eye model structure
%
% Syntax
%  fovea = drasdo.landmarks.fovea( eye )
%
% Description:
%   Calculates the position on the retinal surface of the fovea. Drasdo &
%   Fowler indicate that the fovea is displaced from the longitudinal /
%   optical axis of their model eye, but do not provide the exact values.
%   Watson (2004) provides a numerical approximation to the Drasdo & Fowler
%   result, and assumes that the fovea is displaced 1.5 mm nasally and 0.5
%   mm superiorly on the retinal surface. I implement the Watson fovea in
%   the D&F model eye.

% An initial guess for the fovea is the Watson displacement from the vertex
initialCoords = eye.landmarks.vertex.coords - [0, 1.5, 0.5];

% Obtain the retinal surface quadric
S = eye.retina.S;

% Find the closest point on the quadric surface to the initialCoords
[~,X] = quadric.distancePointEllipsoid( initialCoords', S );
fovea.coords = X';

% Find the geodetic coordinates
fovea.geodetic = quadric.cartToEllipsoidalGeo( fovea.coords', S )';

% Determine where this point falls in visual space by finding the nodal ray
% to this location
[~,fieldAngularPosition] = calcNodalRayToRetina(eye,fovea.coords);

% Save this value
fovea.degField = fieldAngularPosition;


end
