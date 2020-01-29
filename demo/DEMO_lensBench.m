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

% Housekeeping
clear all
close all
clc


%% Plot options
plotOutputRays = false;


%% Set up the refractive index of the medium and lens
mediumRefractiveIndex = 1.0;
lensRefractiveIndex = 2.0;


%% Set the properties of the lenses and position of DLP chip
DLPpostion = 200; % DLP distance from the cornea in mm
diopters = 20; % Lens power
radius = 30; % Lens radius in mm
% Positions of the lenses, in mm.  Must be ordered near to far
lensCenters = [50, 150];    

% Grind the lens
[thickness, curvature] = grindPlus(diopters, radius, lensRefractiveIndex, mediumRefractiveIndex);


%% Initialize the optical system with an eye
% Create an initial optical system in the eyeToCamera (left to right)
% direction with an emmetropic right eye focused at 1.5 meters, with the
% refractive indices for the visible spectrum.
sceneGeometry = createSceneGeometry('spectralDomain','vis');
opticalSystem = sceneGeometry.refraction.retinaToCamera.opticalSystem;

% Insert an iris stop into the system
stopRadius = 2;
rowInsertAfter = find(strcmp(sceneGeometry.refraction.retinaToCamera.surfaceLabels,'lens.front'));
opticalSystem = addIrisStop(opticalSystem,sceneGeometry.eye.stop.center,stopRadius,rowInsertAfter);


%% Add the lenses to the optical system at desired locations
for ll = 1:length(lensCenters)
    opticalSystem = assembleLensWrapper(opticalSystem, lensRefractiveIndex, mediumRefractiveIndex, thickness, curvature, lensCenters(ll), radius);
end


%% Reverse the optical system direction
% The optical system is assembled for ray tracing from left-to-right (away
% from the eye), but we want to do ray tracing from right-to-left (towards
% the eye)
opticalSystem = reverseSystemDirection(opticalSystem);


%% Display the optical system
plotOpticalSystem('surfaceSet',opticalSystem,'addLighting',true,'viewAngle',[0 90]);


%% Add some rays
% Send rays from the horizontal bounds of the DLP chip, which will be 11.5
% mm on either side of the optical axis, and have an angle w.r.t. the
% optical axis of +- 12 degrees
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
        
        % If the outputRay is nan (that is, the ray missed the eye), then
        % extend the ray as it left the last lens surface
        if any(any(isnan(outputRay))) && plotOutputRays
            surfaces = find(~any(isnan(rayPath)));            
            outputRay = rayTraceQuadrics(inputRay, opticalSystem(surfaces,:));
            outputRay(:,2) = outputRay(:,2) * 5; % Pump up the volume
            plotOpticalSystem('newFigure',false,'outputRay',outputRay,'outputRayColor',colors{hh,aa},'viewAngle',[0 90]);            
        end
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
myObj = @(p) norm(diopters-calcDiopters(mySystem(p),true,[-100 -100],radius/3));

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


function z = radiusAtX(F,x)
    myObj = @(z) F(x,0,z);
    z = abs(fzero(myObj,0));
end


function opticalSystemOut = addIrisStop(opticalSystemIn, stopCenter, stopRadius, rowInsertAfter)

S = quadric.scale(quadric.unitSphere,[1,20,20]);
t = stopCenter; t(1) = t(1)-1;
S = quadric.translate(S,t);

stopFront = stopCenter(1);

% Find the x-axis position at which the height of the ellipsoid is equal to
% the desired stop radius
F = quadric.vecToFunc(S);
myObj = @(x) radiusAtX(F,stopFront-1+x)-stopRadius;
x = fzero(myObj,0.9999);

% Create the bounding box
bb = [stopFront-x stopFront  -stopRadius stopRadius -stopRadius stopRadius];

% Assemble a line for the optical system
stopLine = nan(1,19);
stopLine(1:10) = quadric.matrixToVec(S);
stopLine(11) = 1; % Rays traveling left-to-right encounter concave surface
stopLine(12:17) = bb;
stopLine(18) = 1; % Must intersect
% Use the same index of refraction as the surrounding medium, so that the
% aperture stop does not have any effect of refraction.
stopLine(19) = opticalSystemIn(rowInsertAfter,19);

% Add this line to the opticalSystem
opticalSystemOut = [ ...
    opticalSystemIn(1:rowInsertAfter,:); ...
    stopLine; ...
    opticalSystemIn(rowInsertAfter+1:end,:)];

end


