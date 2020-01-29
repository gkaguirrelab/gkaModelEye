function p = plotSurface(S,boundingBox,surfColor,surfAlpha,betaLineColor,omegaLineColor,lineAlpha,bbTol,meshGridSamples)
% Add a 3D plot of the quadric surface to the active figure
%
% Syntax:
%  p = quadric.plotSurface(S,boundingBox,surfColor,surfAlpha,betaLineColor,omegaLineColor,lineAlpha,bbTol)
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
%                           of the surface from 0 (opaque) to 1
%                           (fully transparent)
%   betaLineColor         - 1x3 vector or char string that specifies a
%                           valid MATLAB plor color. E.g. [0.5 0.5 0.5] or
%                           'red'. Used to plot the ellipsoidal geodetic
%                           lines of constant beta. If not defined, these
%                           lines are not shown.
%   omegaLineColor        - 1x3 vector or char string that specifies a
%                           valid MATLAB plor color. E.g. [0.5 0.5 0.5] or
%                           'red'. Used to plot the ellipsoidal geodetic
%                           lines of constant omega. If not defined, these
%                           lines are not shown.
%   bbTol                 - Scalar. Defines the tolerance within which the
%                           intersection must be within the boundingBox.
%                           Default value is 0.1. Handles the situation in
%                           which the intersection is right on the boundary
%                           but is numerically outside.
%   meshGridSamples       - Scalar. Defines the density of surface samples.
%
% Outputs:
%   p                     - Handle to the surface plot object
%
% Examples:
%{
    % Show an ellipsoidal surface with lines of constant beta and omega
    S = quadric.scale(quadric.unitSphere,[4 5 3]);
    boundingBox = [-50 50 -30 30 -20 20];
    figure
    quadric.plotSurface(S, boundingBox, 'k', 0.25, 'red');
    camlight
%}

% Handle incomplete input arguments
if nargin==1
    % synthesize a bounding box that holds the radii of the quadric
    center = quadric.center(S);
    radii = quadric.radii(S);
    boundingBox = [...
        center(1)-radii(1), center(1)+radii(1), ...
        center(2)-radii(2), center(2)+radii(2), ...
        center(3)-radii(3), center(3)+radii(3)];        
    surfColor=[0.9 0.9 0.9];
    surfAlpha=0.8;
    betaLineColor = [];
    omegaLineColor = [];
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==2
    surfColor=[0.9 0.9 0.9];
    surfAlpha=0.8;
    betaLineColor = [];
    omegaLineColor = [];
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==3
    surfAlpha=0.8;
    betaLineColor = [];
    omegaLineColor = [];
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==4
    betaLineColor = [];
    omegaLineColor = [];
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==5
    omegaLineColor = [];
    lineAlpha = 1;
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==6
    lineAlpha = 1;
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==7
    bbTol = 1e-2;
    meshGridSamples = 100;
end

if nargin==8
    meshGridSamples = 100;
end

% Create a linear meshgrid within the boundingBox range
[xx, yy, zz]=meshgrid( linspace(boundingBox(1)-bbTol,boundingBox(2)+bbTol,meshGridSamples),...
    linspace(boundingBox(3)-bbTol,boundingBox(4)+bbTol,meshGridSamples),...
    linspace(boundingBox(5)-bbTol,boundingBox(6)+bbTol,meshGridSamples));

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

% Store the curent hold state
holdState = ishold;

% Plot lines of constant beta
if ~isempty(betaLineColor)
    hold on
    for beta = -90:10:90
        coords =[];
        for omega = -180:3:180
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        inBounds = logical( ...
            (coords(:,1) > boundingBox(1)) .* ...
            (coords(:,1) < boundingBox(2)) .* ...
            (coords(:,2) > boundingBox(3)) .* ...
            (coords(:,2) < boundingBox(4)) .* ...
            (coords(:,3) > boundingBox(5)) .* ...
            (coords(:,3) < boundingBox(6)));
        if ~isempty(coords(inBounds,:))
            lh = plot3(coords(inBounds,1),coords(inBounds,2),coords(inBounds,3),'LineStyle','-','Color',betaLineColor);
            lh.Color(4) = lineAlpha;
        end
    end
    if ~holdState
        hold off
    end
end

% Plot lines of constant omega
if ~isempty(omegaLineColor)
    hold on
    for omega = -180:10:180
        coords =[];
        for beta = -90:3:90
            coords(end+1,:)=quadric.ellipsoidalGeoToCart( [beta, omega, 0], S );
        end
        inBounds = logical( ...
            (coords(:,1) > boundingBox(1)) .* ...
            (coords(:,1) < boundingBox(2)) .* ...
            (coords(:,2) > boundingBox(3)) .* ...
            (coords(:,2) < boundingBox(4)) .* ...
            (coords(:,3) > boundingBox(5)) .* ...
            (coords(:,3) < boundingBox(6)));
        if ~isempty(coords(inBounds,:))
            lh = plot3(coords(inBounds,1),coords(inBounds,2),coords(inBounds,3),'LineStyle','-','Color',omegaLineColor);
            lh.Color(4) = lineAlpha;
        end
    end
    if ~holdState
        hold off
    end
end

end