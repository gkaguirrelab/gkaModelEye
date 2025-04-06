%% t12_eyePoseSearch
% Illustrate a search across eyePose values to fit eye features.
%{
11: Tweet text goes here
%}
%


% Create some default eye features
eyePose = [8 -5 0 2.5];
cameraTrans = [-4; -3; 0];
sceneGeometry=createSceneGeometry();
[ targetEllipse, glintCoord ] = projectModelEye(eyePose,sceneGeometry,'cameraTrans',cameraTrans);
[ Xp, Yp ] = ellipsePerimeterPoints( targetEllipse, 10 );

% Save the true values to display on the gif
xTrue = [eyePose, cameraTrans'];

% Define some inaccurate x0 values, to throw the search initialization off
% the trail for a better demo
cameraTransX0 = [0; 0; 0];
eyePoseX0 = [0 0 0 1];


% Search and save the search history
[~, ~, ~, ~, ~, ~, xHist] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, 'eyePoseX0', eyePoseX0, 'cameraTransX0', cameraTransX0);

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = 'demo/t12_eyePoseSearch.gif';

% These are the elements of the model eye that we will render
modelEyeLabelNames = {'retina' 'pupilEllipse' 'cornea' 'glint_01'};
modelEyePlotColors = {'.w' '-g' '.y' '*r'};

% Get the image dimensions
imDims = sceneGeometry.cameraIntrinsic.sensorResolution;

% Initialize the gif with frames that show the features to be fit
renderEyePose([0 0 0 1],sceneGeometry,'newFigure', true,'modelEyeLabelNames',{});
hold on
gif(gifSaveName);
text(10,450,sprintf('eyePose_true = [ % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f ]',xTrue), 'Color', 'g','Interpreter','none');
for ii = 1:5; gif; end
t1 = text(50,50,'Eye features to be fit','FontSize',20);
for ii = 1:30; gif; end
scatter(Xp,Yp,25,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
t2 = text(50,100,'pupil perimeter','FontSize',14);
for ii = 1:30; gif; end
delete(t2);
for ii = 1:5; gif; end
scatter(glintCoord(1),glintCoord(2),75,'o','MarkerEdgeColor', 'r');
t3 = text(50,150,'glint','FontSize',14);
for ii = 1:30; gif; end
delete(t1);
delete(t3);
for ii = 1:30; gif; end


%% Loop over the search history
plotHandles = [];
poseText = [];
for ii = 1:size(xHist,1)

    if ii ~=size(xHist,1) && ii>1 && max(abs(xHist(ii,:)-xHist(ii-1,:)))<0.001
        continue
    end
    
    delete(plotHandles)
    delete(poseText)
        [~, plotHandles] = renderEyePose(xHist(ii,1:4), sceneGeometry, 'newFigure', false, ...
            'cameraTrans',xHist(ii,5:7)',...
            'modelEyeAlpha',0.5,...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
    
    poseText = text(10,420,sprintf('eyePose__fit  = [ % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f ]',xHist(ii,:)),'Interpreter','none');
        
    % This updates the gif
    gif; gif;
    
end

poseText = text(10,420,sprintf('eyePose__fit  = [ % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f, % 2.1f ]',xHist(ii,:)), 'Color', 'g','Interpreter','none');

for ii = 1:60; gif; end

% Close any open windows
close all

