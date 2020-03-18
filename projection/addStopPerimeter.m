function [eyePoints, pointLabels] = addStopPerimeter(sceneGeometry,p,eyePose)

%% Define an eye in eye coordinates
% This coordinate frame is in mm and has the dimensions (p1,p2,p3). The
% diagram is of a cartoon pupil, viewed directly from the front.
%
% Coordinate [0,0,0] corresponds to the apex (front surface) of the cornea,
% with the model eye having the property of the optical and pupil axes of
% the eye being aligned. The first dimension is depth, and has a negative
% value toward the back of the eye.
%
%                 |
%     ^         __|__
%  +  |        /     \
% p3  -  -----(   +   )-----
%  -  |        \_____/
%     v           |
%                 |
%
%           - <--p2--> +
%
% For the right eye, negative values on the p2 dimension are more temporal,
% and positive values are more nasal. Positive values of p3 are upward,
% and negative values are downward


stopRadius = eyePose(4);

nStopPerimPoints = p.Results.nStopPerimPoints;

%% Define points around the elliptical aperture stop

% The eccentricity of the pupil aperture is given by a stored function
stopEccenFunc = str2func(sceneGeometry.eye.stop.eccenFcnString);

% Determine the parameters of the ellipse that defines the aperture stop in
% the plane of the iris. The absolute value of stopEccenFunc gives the
% eccentricity. The theta of the stop switches from horizontal to vertical
% when the stop passes from a negative to positive eccentricity, passing
% through circular at an eccentricity of 0.
stopEccen = stopEccenFunc(stopRadius);
stopEllipse = [sceneGeometry.eye.stop.center(2) , ...
    sceneGeometry.eye.stop.center(3), ...
    pi*stopRadius^2, ...
    abs(stopEccen),...
    sceneGeometry.eye.stop.thetas(1+(stopEccen>0))];

% Obtain the points on the perimeter of the ellipse
[p2p, p3p] = ellipsePerimeterPoints( stopEllipse, nStopPerimPoints, p.Results.stopPerimPhase );

% Place these points into the eyeWorld coordinates. Optionally create
% separate front and back stop perimeters to model iris thickness.
if sceneGeometry.eye.iris.thickness~=0
    stopPoints = zeros(nStopPerimPoints*2,3);
    stopPoints(1:nStopPerimPoints*2,3) = [p3p; p3p];
    stopPoints(1:nStopPerimPoints*2,2) = [p2p; p2p];
    stopPoints(1:nStopPerimPoints,1) = sceneGeometry.eye.stop.center(1)+sceneGeometry.eye.iris.thickness/2;
    stopPoints(nStopPerimPoints+1:nStopPerimPoints*2,1) = sceneGeometry.eye.stop.center(1)-sceneGeometry.eye.iris.thickness/2;
    
    % Create labels for the stopPerimeter points
    tmpLabelsFront = repmat({'stopPerimeterFront'},nStopPerimPoints, 1);
    tmpLabelsBack = repmat({'stopPerimeterBack'},nStopPerimPoints, 1);
    pointLabels = [tmpLabelsFront; tmpLabelsBack];
else
    stopPoints = zeros(nStopPerimPoints,3);
    stopPoints(1:nStopPerimPoints,3) = p3p;
    stopPoints(1:nStopPerimPoints,2) = p2p;
    stopPoints(1:nStopPerimPoints,1) = sceneGeometry.eye.stop.center(1);
    
    % Create labels for the stopPerimeter points
    tmpLabels = repmat({'stopPerimeter'},nStopPerimPoints, 1);
    pointLabels = tmpLabels;
end
eyePoints = stopPoints;


end