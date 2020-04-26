function opticalSystemOut = addBiconvexLens( opticalSystemIn, lensCenter, diopters, lensRadius, lensRefractiveIndex, mediumRefractiveIndex, surfaceLabel, surfaceColor )
% Add a spheric biconvex (plus) lens to the optical system
%
% Syntax:
%   opticalSystemOut = addBiconvexLens( opticalSystemIn, lensCenter, diopters, radius, lensRefractiveIndex, mediumRefractiveIndex, surfaceLabel, surfaceColor )
%
% Description:
%   Adds a positive diopter, spheric double convex lens to the optical
%   system, centered on the optical axis. The opticalSystem input can be a
%   struct or matrix, but must be in the eyeToCamera direction.
%
%   Note that because this is a simple spheric lens, it will suffer from
%   spherical aberation. Therefore, the assignment of power to this lens is
%   only approximate.
%
% Inputs:
%   opticalSystemIn       - Struct or matrix. If struct, must have the
%                           fields {opticalSystem, surfaceLabels,
%                           surfaceColors}. The matrix form is mx19. See
%                           "assembleOpticalSystem.m" for details. 
%   lensCenter            - Scalar. The position of the stop on the optical 
%                           axis in mm. Defaults to -3.9 mm.
%   diopters              - Scalar. The power of the lens.
%   lensRadius            - Scalar. Radius of the lens in mm.
%   lensRefractiveIndex   - Scalar. Defaults to 2.
%   mediumRefractiveIndex - Scalar. Defaults to 1.
%   surfaceLabel          - Char vector or string. The label for the lens.
%   surfaceColor          - 3x1 vector that specifies the color to be used
%                           to display the lens in the plotOpticalSystem
%                           function.
%
% Outputs:
%   opticalSystemOut      - Struct or matrix. Will match the type of
%                           opticalSystemIn
%
% Examples:
%{
%}


%% Handle nargin
switch nargin
    case 1:3
        error('addBiconvexLens:inputError','Insufficient input arguments');
    case 4
        lensRefractiveIndex = 2;
        mediumRefractiveIndex = 1;
        surfaceLabel = 'biconvexLens';
        surfaceColor = [0.5 0.5 0.5];
    case 5
        mediumRefractiveIndex = 1;
        surfaceLabel = 'biconvexLens';
        surfaceColor = [0.5 0.5 0.5];
    case 6
        surfaceLabel = 'biconvexLens';
        surfaceColor = [0.5 0.5 0.5];
    case 7
        surfaceColor = [0.5 0.5 0.5];
    case 8
        % All is well
    otherwise
        error('addBiconvexLens:inputError','Too many input arguments');
end


%% Obtain the radius and curvature for the specified lens power
[thickness, curvature] = grindPlus(diopters, lensRadius, lensRefractiveIndex, mediumRefractiveIndex);


%% Prepare the optical system

% Extract the opticalSystem matrix
opticalSystemOut = opticalSystemIn;
if isstruct(opticalSystemIn)
    opticalSystemMatrix = opticalSystemIn.opticalSystem;
else
    opticalSystemMatrix = opticalSystemIn;
end

% Test that we have a matrix with a valid system direction
if ~isempty(opticalSystemMatrix) && ~contains(calcSystemDirection(opticalSystemMatrix),'eyeToCamera')
    error('addBiconvexLens:invalidSystemDirection','Lenses are only added to an optical system in the eyeToCamera direction')
end

% Strip the optical system of nan rows. We will add these back before
% returning
numRows = size(opticalSystemMatrix,1);
opticalSystemMatrix = opticalSystemMatrix(sum(isnan(opticalSystemMatrix),2)~=size(opticalSystemMatrix,2),:);


%% Add the lens
backCenter = lensCenter+(thickness/2) + (-curvature);
frontCenter = lensCenter-(thickness/2) + (curvature);
opticalSystemMatrix = assembleLensSystem(opticalSystemMatrix, lensRefractiveIndex, mediumRefractiveIndex, -curvature, backCenter, curvature, frontCenter, lensCenter, lensRadius);


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
if numRows > size(opticalSystemMatrix,1)
    opticalSystemMatrix = [opticalSystemMatrix; ...
        nan(numRows-size(opticalSystemMatrix,1),19)];
end


%% Handle labels, colors, and the output variable
if isstruct(opticalSystemIn)
    opticalSystemOut.opticalSystem = opticalSystemMatrix;
    
    surfaceLabels = opticalSystemIn.surfaceLabels;
    surfaceLabels = [surfaceLabels; ...
        [surfaceLabel '.back']; ...
        [surfaceLabel '.front']];
    surfaceColors = opticalSystemIn.surfaceColors;
    surfaceColors = [surfaceColors; ...
        surfaceColor; ...
        surfaceColor];
    opticalSystemOut.surfaceLabels = surfaceLabels;
    opticalSystemOut.surfaceColors = surfaceColors;
else
    opticalSystemOut = opticalSystemMatrix;
end

end



%% LOCAL FUNCTIONS

function opticalSystemOut = assembleLensWrapper(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, thickness, curvature, lensCenter, radius)
% Adds a positive diopter, double convex lens to the optical system in the
% eyeToCamera direction. This wrapper function allows the lens to be
% specified in terms of thickness, radius, and lens center

backCenter = lensCenter+(thickness/2) + (-curvature);
frontCenter = lensCenter-(thickness/2) + (curvature);

opticalSystemOut = assembleLensSystem(opticalSystemIn, lensRefractiveIndex, mediumRefractiveIndex, -curvature, backCenter, curvature, frontCenter, lensCenter, radius);
end


function [thickness, curvature] = grindPlus(diopters, radius, lensRefractiveIndex, mediumRefractiveIndex)
% Implements a non-linear search to create a plus lens with the specified
% properties

% Anonymous function that returns an optical system with the candidate lens
opticalSystem = [nan(1,18) mediumRefractiveIndex];
mySystem = @(p) assembleLensWrapper(opticalSystem, lensRefractiveIndex, mediumRefractiveIndex, p(1), p(2), 0, radius);

% A non-linear constraint the guides the lens to have the desired radius
myConstraint = @(p) checkLensShape(mySystem(p),radius);

% The objective, which is the requested lens power in diopters
% NOTE THE FOLLOWING HACK -- We are currently creating bi-convex, spherical
% lenses. These are subject to spherical aberration. As a consequence, the
% power of the lens varies depending upon the distance from the optical
% axis at which the power is evaluated. The objective function evaluates
% the lens power at a height of 0.65*radius, as this is observed to produce
% lenses with a central tendency of the focal length that matches the
% called-for power. Really, I should be designing aspheric lenses here, but
% haven't yet implemented this.
myObj = @(p) norm(diopters-calcOpticalPower(mySystem(p),true,[-100 -100],radius*0.65));

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