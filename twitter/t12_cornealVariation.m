%% t12_cornealVariation
% Illustrate the variation in the shape of the cornea
%{
12: The shape of the cornea is well approximated by an ellipsoid. There is individual variation in the ellipsoid parameters across subjects.
%}
%
%


% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = '~/Desktop/t12_cornealVariation.gif';

% The range of cornea parameters
k(1,:) = linspace(7.1450,18.3152,5);
k(2,:) = linspace(40.8458,43.1315,5);
k(3,:) = linspace(48.3340,18.3152,5);
k(4,:) = [45 22.5 0 22.5 45];
k(5,:) = linspace(1.8247,2.5619,5);
k(6,:) = linspace(-0.9422,0.8724,5);


% Plot the initial state
corneaAxialRadius = k(1,3);
kvals = k(2:6,3);
sceneGeometry = createSceneGeometry('kvals',kvals','corneaAxialRadius',corneaAxialRadius);
[figHandle, plotHandles] = plotOpticalSystem('surfaceSet',sceneGeometry.refraction.stopToMedium.opticalSystem([1 3],:),'addLighting',true,'newFigure',true);
plot3([-7 0.5],[0 0],[0 0],'-r');
xlim([-7 0.5])
ylim([-10 10])
zlim([-10 10])
set(gcf,'color','w');
drawnow

% This command opens the gif object
gif(gifSaveName,'frame',gcf);

%% Loop over corneal parameters
for pp = 1:6
    
    
    corneaAxialRadius = k(1,3);
    kvals = k(2:6,3);
    
    for vv = [2 1 2 3 4 5 5 4]
        if pp==1
            corneaAxialRadius = k(1,vv);
        else
            kvals(pp) = k(pp,vv);
        end
        
        % Update the plot
        sceneGeometry = createSceneGeometry('kvals',kvals','corneaAxialRadius',corneaAxialRadius);
        delete(plotHandles)
        [~, plotHandles] = plotOpticalSystem('surfaceSet',sceneGeometry.refraction.stopToMedium.opticalSystem([1 3],:),'addLighting',false,'newFigure',false);
        xlim([-7 0.5])
        ylim([-10 10])
        zlim([-10 10])
        drawnow
        
        % Update the gif
        gif
        
    end
    
end

% Close any open windows
close all

