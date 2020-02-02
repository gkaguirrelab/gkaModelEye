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


%% Set the properties of the lenses and position of DLP chip
DLPpostion = 200; % DLP distance from the cornea in mm
diopters = 20; % Lens power
radius = 30; % Lens radius in mm
% Positions of the lenses, in mm.  Must be ordered near to far
lensCenters = [50, 150];    


%% Initialize the optical system with an eye
% Create an initial optical system in the eyeToCamera (left to right)
% direction with an emmetropic right eye focused at 1.5 meters, with the
% refractive indices for the visible spectrum.
sceneGeometry = createSceneGeometry('spectralDomain','vis');
opticalSystemStruct = sceneGeometry.refraction.retinaToCamera;

% Insert an iris stop into the system
irisApertureRadius = 2;
opticalSystemStruct = addStopAfter(opticalSystemStruct,irisApertureRadius);



%% Add the lenses to the optical system at desired locations
for ll = 1:length(lensCenters)
    opticalSystemStruct = addBiconvexLens( opticalSystemStruct, lensCenters(ll), diopters, radius);
end


%% Reverse the optical system direction
% The optical system is assembled for ray tracing from left-to-right (away
% from the eye), but we want to do ray tracing from right-to-left (towards
% the eye)
opticalSystem = opticalSystemStruct.opticalSystem;
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
        
    end
end




