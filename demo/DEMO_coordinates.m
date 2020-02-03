%% Demo coordinate axes
% Plot and label the axes used in the code and projections
% 
% Description
%   The representation of the eye and the scene in the model code makes use
%   of two different coordinate frames. This script creates a labeled 3D
%   plot that illustrates these frames and their relationship.
%


% Housekeeping
close all
clear
clc


sceneGeometry = createSceneGeometry();

figure


subplot(1,2,1)

% Add the eye
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'newFigure',false,'surfaceAlpha',0.4);

% Add axis lines for the eyeWorld space
plot3([-25 25],[0 0],[0 0],'-k')
plot3([0 0],[-25 25],[0 0],'-k')
plot3([0 0],[0 0],[-25 25],'-k')

% Add end-axis labels
text(-27,0,0,'-p1','HorizontalAlignment','center'); text(+27,0,0,'+p1','HorizontalAlignment','center');
text(0,-27,0,'-p2','HorizontalAlignment','center'); text(0,+27,0,'+p2','HorizontalAlignment','center');
text(0,0,-27,'-p3','HorizontalAlignment','center'); text(0,0,+27,'+p3','HorizontalAlignment','center');

% Indicate floor
f = fill3([-35 35 -35 35],[35 35 -35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)
f = fill3([-35 -35 35 35],[35 -35 35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)

% Handle title and axis
title('eye coordinates')
axis off



subplot(1,2,2)

% Add the eye
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'newFigure',false,'surfaceAlpha',0.4);

% Add axis lines for the eyeWorld space
plot3([-25 25],[0 0],[0 0],'-r')
plot3([0 0],[-25 25],[0 0],'-r')
plot3([0 0],[0 0],[-25 25],'-r')

% Add end-axis labels
text(-27,0,0,'-z','HorizontalAlignment','center'); text(+27,0,0,'+z','HorizontalAlignment','center');
text(0,-27,0,'-x','HorizontalAlignment','center'); text(0,+27,0,'+x','HorizontalAlignment','center');
text(0,0,-27,'-y','HorizontalAlignment','center'); text(0,0,+27,'+y','HorizontalAlignment','center');

% Indicate floor
f = fill3([-35 35 -35 35],[35 35 -35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)
f = fill3([-35 -35 35 35],[35 -35 35 -35],[-35 -35 -35 -35],[0.7 0.4 0.1],'LineStyle','none');
set(f,'facealpha',0.05)

% Handle title and axis
title('world coordinates')
axis off

