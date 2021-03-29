function p = plotImplicitSurface(S,boundingBox,surfColor,surfAlpha,betaLineColor,omegaLineColor,lineAlpha)
% Add a 3D plot of the quadric surface to the active figure
%
% Syntax:
%  p = quadric.plotImplicitSurface(S,boundingBox,surfColor,surfAlpha,betaLineColor,omegaLineColor,lineAlpha)
%
% Description:
%   Plot the quadric surface as an implicit function. This is more precise
%   than plotGridSurface, but the resulting plot object is heavy weight and
%   slow to manipulate.
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
%                           (fully transparent).
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
%   lineAlpha             - Scalar, range 0-1. Specifies the transparency
%                           of the beta and omega lines from 0 (opaque) to
%                           1 (fully transparent).
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
    quadric.plotImplicitSurface(S, boundingBox, 'k', 0.25, 'red');
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
    lineAlpha = 1;
end

if nargin==2
    surfColor=[0.9 0.9 0.9];
    surfAlpha=0.8;
    betaLineColor = [];
    omegaLineColor = [];
    lineAlpha = 1;
end

if nargin==3
    surfAlpha=0.8;
    betaLineColor = [];
    omegaLineColor = [];
    lineAlpha = 1;
end

if nargin==4
    betaLineColor = [];
    omegaLineColor = [];
    lineAlpha = 1;
end

if nargin==5
    omegaLineColor = [];
    lineAlpha = 1;
end

if nargin==6
    lineAlpha = 1;
end


% Plot and define plot appearance
p = fimplicit3(quadric.vecToFunc(S),boundingBox);
p.FaceColor = surfColor;
p.FaceAlpha = min([1 surfAlpha]);
p.EdgeColor = 'none';
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
            if ~strcmp(betaLineColor,'none')
                lh.Color(4) = lineAlpha;
            end
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
            if ~strcmp(omegaLineColor,'none')
            lh.Color(4) = lineAlpha;
            end
        end
    end
    if ~holdState
        hold off
    end
end

end