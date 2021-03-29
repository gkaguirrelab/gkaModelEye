function coordinates = surfaceGrid(S, boundingBox, meshGridSamples, gridType, bbTol)
% Returns a set coordinates evenly distributed on the quadric surface
%
% Syntax:
%  coordinates = quadric.surfaceGrid(S, boundingBox, meshGridSamples, gridType, bbTol)
%
% Description:
%   Returns a set of coordinates that are spaced across the quadric
%   surface. The gridType variable defines how the points are spaced, with
%   the options being: 'parametricPolar', 'ellipsoidalPolar', or
%   'cartesian'.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates are reported.
%   meshGridSamples       - Scalar. The density of the mesh grid.
%   gridType              - Char vector / string. Valid values are:
%                              'cartesian' - the default
%                              'parametricPolar' - non orthogonal lat/long
%                              'ellipsoidalPolar' - Jacobian, orthogonal
%   bbTol                 - Scalar. Defines the tolerance within which the
%                           intersection must be within the boundingBox.
%                           Default value is 0.01. Handles the situation in
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
    S = quadric.scale(quadric.unitSphere,[40 15 30]);
    boundingBox = [0 50 -30 30 -20 20];
    coordinates = quadric.surfaceGrid(S,boundingBox);
    plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'.r')
    axis equal
%}
%{
    eye = modelEyeParameters();
    coordinates = quadric.surfaceGrid(...
        eye.cornea.front.S,...
        eye.cornea.front.boundingBox,...
        23, ...
        'parametricPolar');
    plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'.r')
    axis equal
%}

% Handle incomplete input arguments
if nargin == 2
    meshGridSamples = 100;
    gridType = 'cartesian';
    bbTol = 1e-2;
end

if nargin == 3
    gridType = 'cartesian';
    bbTol = 1e-2;
end

if nargin==4
    bbTol = 1e-2;
end

% Obtain the surfaceGrid
switch gridType

    case 'parametricPolar'

        % Assemble the set of coordinates
        coordinates = [];
        for latitude = linspace(-90,90,meshGridSamples)
            for longitude=linspace(-180,180,meshGridSamples)
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
        for latitude = linspace(-90,90,meshGridSamples)
            for longitude=linspace(-180,180,meshGridSamples)
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

    case 'cartesian'
        
        % Create a linear meshgrid within the boundingBox range
        [xx, yy, zz]=meshgrid( linspace(boundingBox(1)-bbTol,boundingBox(2)+bbTol,meshGridSamples),...
            linspace(boundingBox(3)-bbTol,boundingBox(4)+bbTol,meshGridSamples),...
            linspace(boundingBox(5)-bbTol,boundingBox(6)+bbTol,meshGridSamples));
        
        % Obtain the polynomial function for the quadric surface
        F = quadric.vecToFunc(S);
        
        % Find the vertices that are on the quadric surface
        vf = isosurface(xx, yy, zz, F(xx, yy, zz), 0);
        coordinates = vf.vertices;
        
    otherwise
        error('surfaceGrid:undefinedGridType','Not a valid gridType');
end

end % surfaceGrid

