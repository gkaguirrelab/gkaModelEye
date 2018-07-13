function coordVals = surfaceGrid(S,boundingBox,vertexDensity, polarMesh, bbTol)
%
% Inputs:
%   S                     - A quadric surface in either 1x10 vector form or
%                           4x4 matrix form.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates of an intersection are reported.
%                           Intersection points outside of these bounds are
%                           returned as nans. If not defined no bounds will
%                           be enforced.
%   vertexDensity         - Scalar. The density of the mesh grid. Defaults
%                           to unity.
%   polarMesh             - Logical. If set to true, produces a polar
%                           distribution of points. This is valid only if
%                           the quadric is an ellipsoid.
%
% Outputs:
%   vertices              - nx3 matrix that provides the [x; y; z]
%                           coordinates of n points that are on the surface
%                           of the quadric within the boundingBox
%   faces                 - fx3 matrix that provides the faces of the
%                           surface
%
% Examples:
%{
    S = quadric.unitTwoSheetHyperboloid;
    S = quadric.scale(S,[5 5 1]);
    boundingBox = [0 50 -30 30 -20 20];
    vertices = quadric.surfaceMeshVertices(S,boundingBox,50);
    plot3(vertices(:,1),vertices(:,2),vertices(:,3),'.r')
%}

% Handle incomplete input arguments
if nargin == 2
    vertexDensity = 1;
end

if nargin == 3
    polarMesh = true;
    bbTol = 1e-2;
end

if nargin==4
    bbTol = 1e-2;
end

% If the quadric surface was passed in matrix form, convert to vec
if isequal(size(S),[4 4])
    S = quadric.matrixToVec(S);
end

% Obtain the surfaceGrid
if polarMesh
    % Shift the quadric to the center, and the boundingBox by the
    % corresponding amount
    Xt = quadric.center(S);
    St = quadric.translate(S,-Xt);
    boundingBox(1:2) =  boundingBox(1:2)-Xt(1);
    boundingBox(3:4) =  boundingBox(3:4)-Xt(2);
    boundingBox(5:6) =  boundingBox(5:6)-Xt(3);
    % Obtain the radii.
    radii = quadric.radii(St);
    coordVals = [];
    for latitude = linspace(-180,180,vertexDensity*2)
        for longitude=linspace(-180,180,vertexDensity)
            X = quadric.geodeticToCart( [latitude; longitude; 0], radii );
            % Store the coordinate value.
            coordVals(end+1,:)= X;
        end
    end
    % Re-order the coordVal dimensions to match the axis order for the
    % quadric.
    axisOrder = quadric.axisOrder(S);
    coordVals = coordVals(:,axisOrder);    
    % Remove the coordinates that are outside the bounding box
    retainCoords = (coordVals(:,1) > boundingBox(1)-bbTol) .* (coordVals(:,1) < boundingBox(2)+bbTol) .* ...
        (coordVals(:,2) > boundingBox(3)-bbTol) .* (coordVals(:,2) < boundingBox(4)+bbTol) .* ...
        (coordVals(:,3) > boundingBox(5)-bbTol) .* (coordVals(:,3) < boundingBox(6)+bbTol);
    coordVals = coordVals(logical(retainCoords),:);
    % Shift the coord vals to be w.r.t. the original quadric center
    coordVals = coordVals + Xt';
else
    coordVals = [linspace(boundingBox(1),boundingBox(2),vertexDensity); ...
	linspace(boundingBox(3),boundingBox(4),vertexDensity); ...
    	linspace(boundingBox(5),boundingBox(6),vertexDensity)];
end

end

