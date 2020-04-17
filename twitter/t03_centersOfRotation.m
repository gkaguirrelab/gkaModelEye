%% t03_centersOfRotation
% Show on a schematic illustration the center of rotation for horizontal
% and vertical eye movements
%{
03: The eye rotates. Fry & Hill (1962, 1963) found that the center of rotation differs by rotation direction, and is slightly deeper (and shifted nasally from the optical axis) for horizontal eye movements (*) than vertical eye movements (o).
%}

% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another
sceneGeometry=createSceneGeometry();

% In the Connectome dataset we find that the rotation centers are about 10%
% less deep than Fry & Hill reported
%
%       Fry, G. A., and W. W. Hill. "The center of rotation of the
%       eye." Optometry and Vision Science 39.11 (1962): 581-595.
%
%       Fry, Glenn A., and W. W. Hill. "The mechanics of elevating the
%       eye." Optometry and Vision Science 40.12 (1963): 707-716.
%
sceneGeometry.eye.rotationCenters.azi = sceneGeometry.eye.rotationCenters.azi*0.9;
sceneGeometry.eye.rotationCenters.ele = sceneGeometry.eye.rotationCenters.ele*0.9;

% Plot the schematic eye in red, including the rotation centers
plotModelEyeSchematic(sceneGeometry.eye,'view','horizontal','plotRotationCenters',true,'plotColor','r')

% Now over-plot in black, so that the rotation centers stand out
plotModelEyeSchematic(sceneGeometry.eye,'view','horizontal','newFigure',false,'plotColor','k')

% Clean up the plot limits
xlim([-25 5])
ylim([-15 15])
axis square