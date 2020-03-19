function [eyePoints, pointLabels] = addStopPerimeter(sceneGeometry,p,eyePose)
% Initialize the eyePoints with the perimeter of the iris aperture stop
%
% Syntax:
%  [eyePoints, pointLabels] = addStopPerimeter(sceneGeometry,p,eyePose)
%
% Description:
%   Define a set of points in eyeWorld coordinates that describe the
%   eliptical perimeter of the aperture stop of the iris.
%
%	The eyeWorld coordinate frame is in mm and has the dimensions
%	(p1,p2,p3). The diagram is of a cartoon pupil, viewed directly from the
%	front (i.e., from the perspective of a camera looking at the eye of a
%	subject).
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
%   Coordinate [0,0,0] corresponds to the apex (front surface) of the
%   cornea (when the apex of the cornea is aligned with the optical axis).
%   The first dimension is depth, and has a negative value toward the back
%   of the eye. For the right eye, negative values on the p2 dimension are
%   more temporal, and positive values are more nasal. Positive values of
%   p3 are upward, and negative values are downward
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   p                     - Structure. The structure returned by the
%                           parameter parser in the calling function.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%
% Outputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%



% Extract some values for clarity in the code that follows
stopRadius = eyePose(4);
nStopPerimPoints = p.Results.nStopPerimPoints;

% The eccentricity of the aperture stop is given by a stored function. The
% absolute value of this function gives the eccentricity, and the sign is
% used to set theta (the tilt of the ellipse) below.
stopEccenFunc = str2func(sceneGeometry.eye.stop.eccenFcnString);
stopEccen = stopEccenFunc(stopRadius);

% Define here the paramters of the ellipse (in transparent form) that
% describes the aperture stop. These are:
%
% -	The x and y centers of the ellipse (in the p2 and p3 coordinates of the
%   eyeWorld system)
% - The area of the ellipse, which is defined as equal in area to a circle
%   of the specified stopRadius.
% - The non-linear eccentricity of the ellipse
% - The theta (tilt) of the ellipse, which switches from horizontal to
%   vertical when the stop passes from a negative to positive eccentricity,
%   passing through circular at an eccentricity of 0.
%
stopEllipse = [sceneGeometry.eye.stop.center(2) , ...
    sceneGeometry.eye.stop.center(3), ...
    pi*stopRadius^2, ...
    abs(stopEccen),...
    sceneGeometry.eye.stop.thetas(1+(stopEccen>0))];

% Obtain the points on the perimeter of this ellipse
[p2p, p3p] = ellipsePerimeterPoints( stopEllipse, nStopPerimPoints, p.Results.stopPerimPhase );

% Place these points into the eyeWorld coordinates at the depth (p1) of the
% iris plane. Optionally create separate front and back stop perimeters to
% model iris thickness.
if sceneGeometry.eye.iris.thickness~=0
    
    % The aperture stop has a thickness. Create two sets of points at the
    % front and back of the iris plane.
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
    
    % The aperture stop has zero thickness. One set of points.
    stopPoints = zeros(nStopPerimPoints,3);
    stopPoints(1:nStopPerimPoints,3) = p3p;
    stopPoints(1:nStopPerimPoints,2) = p2p;
    stopPoints(1:nStopPerimPoints,1) = sceneGeometry.eye.stop.center(1);
    
    % Create labels for the stopPerimeter points
    tmpLabels = repmat({'stopPerimeter'},nStopPerimPoints, 1);
    pointLabels = tmpLabels;
end

% Store the stop points in the eyePoints variable to return
eyePoints = stopPoints;


end % addStopPerimeter