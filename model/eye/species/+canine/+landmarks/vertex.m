function vertex = vertex( eye )
% Returns the vertex for the landmarks sub-field of an eye model structure
%
% Syntax
%  vertex = human.landmarks.vertex( eye )
%
% Description:
%   Calculates the position on the retinal surface of the posterior segment
%   vertex.
%
% Inputs
%   eye                   - Structure.
%
% Outputs
%   vertex                - Structure with the subfields degField,
%                           geodetic, and coords
%


% Obtain the quadric form of the retinal surface
S = eye.retina.S;

% Eye landmarks are specified as rotations (in degrees) within the eye
% world coordinate frame for azimuth, elevation, and rotation. Axes are
% defined relative to the optical axis, which itself is set to be aligned
% with the p1 dimension of the eye world coordinate frame.
vertex.degField = [0 0 0];
vertex.geodetic = [-90 -90 0];
vertex.coords = quadric.ellipsoidalGeoToCart(vertex.geodetic,S)';



end
