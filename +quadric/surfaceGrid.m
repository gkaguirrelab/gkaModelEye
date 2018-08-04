function coordinates = surfaceGrid(S, boundingBox, vertexDensity, gridType, bbTol)
% Returns a set coordinates evenly distributed on the quadric surface
%
% Syntax:
%  coordinates = qaudric.surfaceGrid(S, boundingBox, vertexDensity, polarGrid, bbTol)
%
% Description:
%   Returns a set of coordinates that are spaced across the quadric
%   surface. If a polarGrid is requested, then the points are evenly spaced
%   in geodetic coordinates. Otherwise, the points are evenly spaced in
%   Cartesian (linear) coordinates.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates are reported.
%   vertexDensity         - Scalar. The density of the mesh grid. Defaults
%                           to unity.
%   gridType              - Char vector / string. Valid values are:
%                              'linear' - the default
%                              'parametricPolar' - non orthogonal lat/long
%                              'ellipsoidalPolar' - Jacobian, orthogonal
%   bbTol                 - Scalar. Defines the tolerance within which the
%                           intersection must be within the boundingBox.
%                           Default value is 0.1. Handles the situation in
%                           which the intersection is right on the boundary
%                           but is numerically outside.
%
% Outputs:
%   coordinates           - nx3 matrix that provides the [x; y; z]
%                           coordinates of n points that are on the surface
%                           of the quadric within the boundingBox
%
% Examples:
%{
    S = quadric.scale(quadric.unitSphere,[4 5 3]);
    boundingBox = [0 50 -30 30 -20 20];
    coordinates = quadric.surfaceGrid(S,boundingBox,200);
    plot3(coordinates(1,:),coordinates(2,:),coordinates(3,:),'.r')
    axis equal
%}

% Handle incomplete input arguments
if nargin == 2
    vertexDensity = 1;
end

if nargin == 3
    bbTol = 1e-2;
    gridType = 'linear';
end

if nargin==4
    bbTol = 1e-2;
end

% Obtain the surfaceGrid
switch gridType
    case 'parametricPolar'
        % Assemble the set of coordinates
        coordinates = [];
        for latitude = linspace(-90,90,vertexDensity)
            for longitude=linspace(-180,180,vertexDensity)
                X = quadric.parametricGeoToCart( [latitude; longitude; 0], S );
                % Store the coordinate value.
                coordinates(end+1,:)= X;
            end
        end
        % Remove the coordinates that are outside the bounding box
        retainCoords = (coordinates(:,1) > boundingBox(1)-bbTol) .* (coordinates(:,1) < boundingBox(2)+bbTol) .* ...
            (coordinates(:,2) > boundingBox(3)-bbTol) .* (coordinates(:,2) < boundingBox(4)+bbTol) .* ...
            (coordinates(:,3) > boundingBox(5)-bbTol) .* (coordinates(:,3) < boundingBox(6)+bbTol);
        coordinates = coordinates(logical(retainCoords),:);
    case 'ellipsoidalPolar'
        % Assemble the set of coordinates
        coordinates = [];
        for latitude = linspace(-90,90,vertexDensity)
            for longitude=linspace(-180,180,vertexDensity)
                X = quadric.ellipsoidalGeoToCart( [latitude; longitude; 0], S );
                % Store the coordinate value.
                coordinates(end+1,:)= X;
            end
        end
        % Remove the coordinates that are outside the bounding box
        retainCoords = (coordinates(:,1) > boundingBox(1)-bbTol) .* (coordinates(:,1) < boundingBox(2)+bbTol) .* ...
            (coordinates(:,2) > boundingBox(3)-bbTol) .* (coordinates(:,2) < boundingBox(4)+bbTol) .* ...
            (coordinates(:,3) > boundingBox(5)-bbTol) .* (coordinates(:,3) < boundingBox(6)+bbTol);
        coordinates = coordinates(logical(retainCoords),:);
    case 'linear'
        % Produce a set of coordinates that are linearly spaced across the
        % boundingBox.
        [X,Y,Z] = meshgrid(linspace(boundingBox(1),boundingBox(2),vertexDensity), ...
            linspace(boundingBox(3),boundingBox(4),vertexDensity), ...
            linspace(boundingBox(5),boundingBox(6),vertexDensity));
        F = quadric.vecToFunc(S);
        [x,y,z] = ind2sub(size(X),find(F(X,Y,Z) < 1e-6));
        coordinates = [x y z]';
    otherwise
        error('Not a valid surface grid type');
end

end % surfaceGrid

