%% Lens Bench Demo
% A simple optical bench explorer
%
% Description:
%   The origin of the optical system [0,0,0] corresponds to the apex of the
%   cornea of the eye. The optical axis of the eye runs along the x
%   dimension. The eye is positioned at negative x values. This routine
%   allows one to specify the position of bi-convex lenses centered along
%   the optical axis. You must add the lenses in the left to right
%   direction. That is, each subsequent lens must have a larger (more
%   positive) position.
%

clear all
close all
clc

% Set up the properties of the medium and the lenses
mediumRefractiveIndex = 1.0;
lensRefractiveIndex = 2;

% Create an initial optical system in the eyeToCamera (left to right)
% direction with an emmetropic right eye focused at 1.5 meters.
sceneGeometry = createSceneGeometry();
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;

% Design a +20 lens
diopters = 20;
radius = 30; % mm
[thickness, curvature] = grindPlus(diopters, radius, lensRefractiveIndex, mediumRefractiveIndex);

% Add it to the optical system at two locations (must be ordered in the eye
% to camera direction)
lensCenter = 50;
opticalSystem = assembleLensWrapper(opticalSystem, lensRefractiveIndex, mediumRefractiveIndex, thickness, curvature, lensCenter, radius);

lensCenter = 150;
opticalSystem = assembleLensWrapper(opticalSystem, lensRefractiveIndex, mediumRefractiveIndex, thickness, curvature, lensCenter, radius);

% Reverse the optical system to be in cameraToEye orientation
opticalSystem = reverseSystemDirection(opticalSystem);

% Display the optical system and ray
plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'viewAngle',[0 90]);

% Send rays from the horizontal bounds of the DLP chip, which will be 11.5
% mm on either side of the optical axis, and have an angle w.r.t. the
% optical axis of +- 12 degrees
DLPpostion = 200;
angles = [-12,0,12];
horizPos = [-11.5,0,11.5];
colors = {'red','red','red';...
    'green','green','green';...
    'blue','blue','blue'};

% Loop over the desired rays and display
for hh = 1:length(horizPos)
    for aa = 1:length(angles)
        % Create the ray
        inputRay = quadric.normalizeRay(quadric.anglesToRay([DLPpostion;horizPos(hh);0],-180+angles(aa),0));

        % Trace it
        [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem);
        
        % Add it to the plot
        plotOpticalSystem('newFigure',false,'rayPath',rayPath,'rayColor',colors{hh,aa},'viewAngle',[0 90]);
    end
end








%% LOCAL FUNCTIONS

function [thickness, curvature] = grindPlus(diopters, radius, lensRefractiveIndex, mediumRefractiveIndex)
% Implements a non-linear search to create a plus lens with the specified
% properties

% Anonymous function that returns an optical system with the candidate lens
opticalSystem = [nan(1,18) mediumRefractiveIndex];
mySystem = @(p) assembleLensWrapper(opticalSystem, lensRefractiveIndex, mediumRefractiveIndex, p(1), p(2), 0, radius);

% A non-linear constraint the guides the lens to have the desired radius
myConstraint = @(p) checkLensShape(mySystem(p),radius);

% The objective, which is the requested lens power in diopters
myObj = @(p) norm(diopters-calcDiopters(mySystem(p),true,[-100 -100]));

% Set an x0 that is a 20 mm thickness and a curvature based on the thin
% lens approximation.
x0 =  [20 abs((mediumRefractiveIndex-lensRefractiveIndex)/(diopters/2)*1000)];

% Set physically possible bounds on the lens
lb = [0 0];
ub = [Inf Inf];

% Silence the search
options = optimset('fmincon');
options.Display='off';

% Perform the search
p = fmincon(myObj,x0,[],[],[],[],lb,ub,myConstraint,options);

% Report the lens properties
thickness = p(1);
curvature = p(2);

end


function [c, ceq] = checkLensShape(opticalSystem,targetRadius)
% Check if the radius of the lens differs from the target radius. I'm sure
% there is an analytic approach to solving for the radius and curvature
% while getting the diopters right, but we are going to brute-force it here
% with a non-linear search.

S = opticalSystem(2,1:10);
F = quadric.vecToFunc(S);
myFun = @(z) F(0,0,z);
radius = abs(fzero(myFun,0));

c = [];
ceq = targetRadius - radius;

end


function opticalSystemOut = assembleLensWrapper(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, thickness, curvature, lensCenter, radius)
% Adds a positive diopter, double convex lens to the optical system in the
% eyeToCamera direction. This wrapper function allows the lens to be
% specified in terms of thickness, radius, and lens center

backCenter = lensCenter+(thickness/2) + (-curvature);
frontCenter = lensCenter-(thickness/2) + (curvature);

opticalSystemOut = assembleLensSystem(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, -curvature, backCenter, curvature, frontCenter, lensCenter, radius);
end


function opticalSystemOut = assembleLensSystem(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, backCurvature, backCenter, frontCurvature, frontCenter, lensCenter, radius)
% Assembles and returns an optical system matrix given input in the
% eyeToCamera direction

% Define a bounding box for the front surface (from left to right)
boundingBoxLens = [frontCenter-frontCurvature lensCenter  -radius radius -radius radius];

% Add the front spectacle surface to the optical system.
SlensFront = quadric.scale(quadric.unitSphere,[frontCurvature frontCurvature frontCurvature]);
SlensFront = quadric.translate(SlensFront,[frontCenter 0 0]);
lensLine = nan(1,19);
lensLine(1:10) = quadric.matrixToVec(SlensFront);
lensLine(11) = -1; % rays intersect convex lens surface
lensLine(12:17) = boundingBoxLens;
lensLine(18) = 1; % must intersect
lensLine(end) = lensRefractiveIndex;
opticalSystemOut = [opticalSystemIn; lensLine];

% Define a bounding box for the back surface
boundingBoxLens = [lensCenter backCenter-backCurvature -radius radius -radius radius];

% Add the back spectacle surface to the optical system.
SlensBack = quadric.scale(quadric.unitSphere,[backCurvature backCurvature backCurvature]);
SlensBack = quadric.translate(SlensBack,[backCenter 0 0]);
lensLine = nan(1,19);
lensLine(1:10) = quadric.matrixToVec(SlensBack);
lensLine(11) = 1; % rays intersect concave lens surface
lensLine(12:17) = boundingBoxLens;
lensLine(18) = 1; % must intersect
lensLine(end) = mediumRefractiveIndex;
opticalSystemOut = [opticalSystemOut; lensLine];

end