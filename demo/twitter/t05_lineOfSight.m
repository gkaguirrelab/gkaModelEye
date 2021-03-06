%% t05_lineOfSight
% Show a schematic of the fovea and line of sight. The identities and names
% of the various axes of the eye is complicated. A good reference is:
%
%   Atchison, D. A. & Smith, G. Optics of the human eye
%   (Butterworth-Heinemann Oxford, 2000).
%{
05: Unlike a camera, optical elements of the eye are not all aligned. The fovea (*) is typically displaced inferotemporally from the retinal apex. The line-of-sight (-) connects the fixation point, center of entrance pupil, and fovea.
%}
% For details see:
%   calc/calcLineOfSightRay.m
%

% Create a sceneGeometry object that describes the eye, a camera, and their
% position relative to one another. Set the key-value to calculate the
% location of the fovea
eye = modelEyeParameters();

% Calculate the line of sight for the eye
rayDestination = eye.landmarks.fovea.coords;
% Find the sight ray to the fovea (i.e., the line of sight axis)
rayPath = calcSightRayToRetina(eye,rayDestination);

% Plot the schematic eye in red
plotModelEyeSchematic(eye,'view','horizontal','plotColor','r', ...
    'rayPath',rayPath);

% Remove the fovea
eye.landmarks=rmfield(sceneGeometry.eye.landmarks,'fovea');

% Now over-plot in black, without the fovea and line of sight
plotModelEyeSchematic(eye,'view','horizontal','newFigure',false,'plotColor','k');

% Clean up the plot limits
xlim([-25 5])
ylim([-15 15])
axis square