function [opticalAxisRay,errors] = calcOpticalAxis(opticalSystem,rayOriginDistance,cameraMedium)
% Returns the optical axis for an opticalSystem
%
% Syntax:
%  [opticalAxisRay,errors] = calcOpticalAxis(opticalSystem,rayOriginDistance,cameraMedium)
%
% Description
%   The optical axis of a system is the ray that enters and exits the
%   optical system along the same, straight line. For a system of simple
%   spherical lenses, this corresponds to the line that passes through the
%   centers of curvature of the surfaces. When these surfaces are centered
%   on the longitudinal axis, then the longitudinal and optical axes are
%   equivalent. As developed by Harris (2010), more complicated optical
%   systems are likely to have a single optical axis, and simple optical
%   systems may have zero or an infinity of such rays:
%
%       Harris, William F. "Optical axes of eyes and other optical
%       systems." Optometry and Vision Science 86.5 (2009): 537-541.
%
%   For an eye with decententered elements, the optical axis diverges
%   somewhat from the longitudinal axis.
%
% Inputs:
%   opticalSystem         - Either an eye structure (from which a
%                           "mediumToRetina" optical system in air will be
%                           derived), or an mx19 matrix, where m is set by
%                           the key value opticalSystemNumRows. Each row
%                           contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from origin of the longitudinal axis.
%                           Assumed to be 1500 mm if not defined.
%
% Outputs:
%   opticalAxisRay        - 3x2 matrix that is the optical axis ray.
%   errors                - 1x1 matrix with the follow error values:
%                             - The L2 norm of the distance between ray,
%                               and the projected position of the output ray
%
% Examples:
%{
    % Find the optical axis of an eye
    eye = modelEyeParameters();
    [opticalAxisRay, errors] = calcOpticalAxis(eye);
%}


arguments
    opticalSystem
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    cameraMedium = 'air'
end


% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);
    end
end

% Strip the optical system of any rows which are all nans
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);

% Obtain the system direction
systemDirection = calcSystemDirection(opticalSystem, rayOriginDistance);

% Define the objective
myObj = @(p) objective(p,opticalSystem,rayOriginDistance);

% Set the p0 and bounds based upon the ray trace direction
switch systemDirection
    case 'cameraToEye'
        p0 = [0,0,-180,-180];
        lb = [-rayOriginDistance*2,-rayOriginDistance*2,-270,-270];
        ub = [rayOriginDistance*2,rayOriginDistance*2,-90,-90];
    case 'eyeToCamera'
        p0 = [0,0,0,0];
        lb = [-rayOriginDistance*2,-rayOriginDistance*2,-90,-90];
        ub = [rayOriginDistance*2,rayOriginDistance*2,90,90];
    otherwise
        error(['Not a valid system direction: ' systemDirection])
end

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search. We are finding the [x,y] ray start point and angles that yield
% the initial ray that is closest to the property of an optical axis
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% Evaluate the objective once more and construct the optical axis to return
[errors, opticalAxisRay] = objective(p,opticalSystem,rayOriginDistance);

end


%% Local function

function [fVal, inputRay] = objective(p,opticalSystem,rayOriginDistance)

% Assemble the ray
inputRay = quadric.normalizeRay(quadric.anglesToRay([rayOriginDistance;p(1);p(2)],p(3),p(4)));

% Perform the ray trace
outputRay = rayTraceQuadrics(inputRay, opticalSystem);

% Project the output ray back to the rayOrigin
X = outputRay(:,1)-norm(inputRay(:,1)-outputRay(:,1)).*outputRay(:,2);

% The L2 norm of the distance between the starting position of the input
% ray, and the projected position of the output ray
fVal = norm(inputRay(:,1)-X);

end

