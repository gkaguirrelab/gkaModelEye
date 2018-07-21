function p = plotSurface(S,boundingBox,surfColor,surfAlpha)
% Add a 3D plot of the quadric surface to the active figure
%
% Syntax:
%  p = quadric.plotSurface(S,boundingBox,surfColor,surfAlpha)
%
% Description:
%   Creates a meshgrid on the quadric surface and plots this as a mesh
%   surface, using the supplied face color and alpha transparency.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates are reported.
%   surfColor             - 1x3 vector or char string that specifies a
%                           valid MATLAB plor color. E.g. [0.5 0.5 0.5] or
%                           'red'.
%   surfAlpha             - Scalar, range 0-1. Specifies the transparency
%                            of the surface from 0 (opaque) to 1
%                            (fully transparent)
%
% Outputs:
%   p                     - Handle to the surface plot object
%

% Define the level of detail of the surface mesh.
meshGridSamples = 100;

% Create a linear meshgrid within the boundingBox range
[xx, yy, zz]=meshgrid( linspace(boundingBox(1),boundingBox(2),meshGridSamples),...
    linspace(boundingBox(3),boundingBox(4),meshGridSamples),...
    linspace(boundingBox(5),boundingBox(6),meshGridSamples));

% Obtain the polynomial function for the quadric surface
F = quadric.vecToFunc(S);

% Find the vertices that are on the quadric surface
vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);

% Plot and define plot appearance
p = patch(vertices);
p.FaceColor = surfColor;
p.EdgeColor = 'none';
alpha(surfAlpha);
daspect([1 1 1])
view(3);
axis tight
axis equal

end