%% t09_rotationalAsymmetry
% Illustrate aspects of the eye that do not have rotational symmetry
%{
09: Modeling torsion matters as (unlike a camera) the eye is not rotationally symmetric. E.g. (weird eye exagerated x4): a) Cornea is generally wider than tall; b) The aperture stop itself is taller than wide when dilated (also with a nasal top tilt) (HJ Wyatt 1995)
%}
% For details see:
%   eye/species/+human/cornea.m
%   eye/species/+human/calcDerivedParams.m
%

% Save location for the GIF. Sorry for the Mac bias.
gifSaveName = 'demo/t09_rotationalAsymmetry.gif';

% How much do we want to exagerate these non-symmetric features?
exagerate = 4;

% The cornea is well modeled as a tri-axial ellipsoid. It is traditional to
% describe the front surface of the cornea in terms of the optical
% refractive power of the surface. Within this standard, the mean corneal
% power (in diopters) is:
meanCornealPower = 44.9356;

% The cornea normally has some astigmatism, usually with power being
% greater along the vertical axis than the horizontal axis. We will call
% this here the k1k2 difference:
k1k2DiffDiopters = 1.3892;

% Now we make a kvals vector that exagerates the k1k2 difference x5
kvals = [meanCornealPower - k1k2DiffDiopters*exagerate, ...
    meanCornealPower + k1k2DiffDiopters*exagerate, ...
    0, 0 0];

% Time to change the eccentricity of the aperture stop. First, we load a
% file called the "derivedParams" that contains parameter values that
% define the behavior of the shape of the stop
filename = fullfile(replace(replace(mfilename('fullpath'),mfilename(),''),'demo/twitter','model'),'eye','species','+human','derivedParams.mat');
load(filename,'derivedParams');

% Exagerate the maximum eccentricity of the stop
derivedParams.stopEccenParams(4) = derivedParams.stopEccenParams(4) * exagerate;

% Create a sceneGeometry with these kvals
sceneGeometry=createSceneGeometry('kvals',kvals,'derivedParams',derivedParams);

% Remove the effect of corneal refraction. To do so, we will hack the
% index of refraction of the aqueous humor, cornea, and tear film to all
% be one.
sceneGeometry.refraction.stopToMedium.opticalSystem(1:3,19) = 1;

% These are the elements of the model eye that we will render. The "pupil"
% is now equivalent to the aperture stop, so we plot it in red.
modelEyeLabelNames = {'retina' 'pupilEllipse' 'cornea'};
modelEyePlotColors = {'.w' '-r' '.y'};

% The angles across which the eye will rotate
rotationValues = [0:0.5:5,5:-0.5:-5,-5:0.5:-0.5];


%% Loop over eyePoses
for ii = 1:length(rotationValues)
    
    eyePose = [rotationValues(ii) 0 0 3];
    
    if ii == 1
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', true, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
        
        % This command opens the gif object
        gif(gifSaveName);
        
    else
        delete(plotHandles(1:end))
        [~, plotHandles] = renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
            'modelEyeLabelNames', modelEyeLabelNames, ...
            'modelEyePlotColors', modelEyePlotColors);
    end
    
    % This updates the gif
    gif
    
end

% Close any open windows
close all
