function [S,fVal] = fitEllipsoidToPoints(pp,s0,sB)
% Fit the radii of a quadric to set of points
%
% Syntax:
%  [s,fVal] = quadric.fitEllipsoidToPoints(pp,s0,sB)
%
% Description:
%   Returns the coordinates of the points of intersection of a ray with a
%   quadric surface.
%
% Inputs:
%   pp                    - 3xn matrix of points to be fit 
%   s0                    - 3x1 vector of initial values for the radii of
%                           the quadric surface.
%   sB                    - 3x1 vector of bounds on the radii of the 
%                           quadric surface.
%
% Outputs:
%   S                     - 1x10 vector of the quadric surface.
%   fVal                  - Scalar. L2 norm of the distance of pp from S
%
% Examples:
%{
    % Find the radii of the image surface for the emmetropic eye at D0
    eye = modelEyeParameters('spectralDomain',555,'accommodation',0);
    principalPoint = calcPrincipalPoint(eye);
    iFP = []; % internal focal points
    a = 30;
    for hh = -a:5:a
        for vv = -a:5:a
            if norm([hh vv])>a
                continue
            end
            iFP(:,end+1) = calcInternalFocalPoint(eye,[hh vv],Inf,eye.landmarks.incidentNode.coords,principalPoint);
        end
    end
    % Get the retina radii and place bounds
    s0 = quadric.radii(eye.retina.S);
    sB = [0; 10; 10];
    % Shift the iFP points to be on the border of the centered retina
    iFP(1,:) = iFP(1,:)-min(iFP(1,:)) - s0(1);
    [SFit,fVal] = quadric.fitEllipsoidToPoints(iFP,s0,sB);
    % Report the difference between s0 and s
    sD = s0 - quadric.radii(SFit);
    fprintf('The difference [h,v] radii between the emmetropic retina and the image surface: [%2.2f, %2.2f]\n',sD(2:3));
%}

% Functions for the quadric surface and the objective
ellipsoid = @(s) quadric.scale(quadric.unitSphere,s);
myObj = @(s) objective(ellipsoid(s),pp);

% Options
options = optimset('fmincon');
options.Display = 'off';

% Search
[s,fVal]=fmincon(myObj,s0,[],[],[],[],s0-sB,s0+sB,[],options);

% Return the quadric surface
S = ellipsoid(s);

end


% Local function for the objective
function fVal = objective(S,pp)

for ii = 1:size(pp,2)
    d(ii) = quadric.distancePointEllipsoid(pp(:,ii),S);
end
fVal = norm(d);

end
