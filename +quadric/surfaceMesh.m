function [vertices, faces] = surfaceMesh(v,boundingBox,vertexDensity)
%
% Inputs:
%   v                     - A quadric surface in either 1x10 vector form.
%                           If the 4x4 matrix form is passed it will be
%                           converted.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates of an intersection are reported.
%                           Intersection points outside of these bounds are
%                           returned as nans. If not defined no bounds will
%                           be enforced.
%   vertexDensity         - Scalar. The density of the mesh grid. Defaults
%                           to unity.
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

% If the quadric surface was passed in matrix form, convert to vec
if isequal(size(v),[4 4])
    v = quadric.matrixToVec(v);
end

% Obtain the meshgrid
[xx, yy, zz]=meshgrid( linspace(boundingBox(1),boundingBox(2),vertexDensity),...
    linspace(boundingBox(3),boundingBox(4),vertexDensity),...
    linspace(boundingBox(5),boundingBox(6),vertexDensity));

% Create a function handle for the quadric surface
F = quadric.vecToFunc(v);

% Obtain the vertices and faces
fv = isosurface(xx, yy, zz, F(xx, yy, zz), 0);
vertices = fv.vertices;
faces = fv.faces;

end

