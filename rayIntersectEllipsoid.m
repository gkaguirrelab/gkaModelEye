function [ coordsOut, curvature_p1p2, curvature_p1p3 ] = rayIntersectEllipsoid( coordsIn, theta, phi, ellipsoidRadii, ellipsoidCenter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ellipsoidRadii = [3 1 2];
ellipsoidCenter = [0 0 0];
rayP0 = [0 0 0];
theta = deg2rad(0);
phi = deg2rad(0);

% Store the sign of the radius values. They radii must have the same sign.
if prod(sign(ellipsoidRadii))==-1
    error('rayTraceCenteredSurfaces:incompatibleConvexity','The radii of the elliptical lens surface must have the same sign.');
end
radiiSign = sign(ellipsoidRadii(1));
% Convert the radii to their absolute values
ellipsoidRadii = abs(ellipsoidRadii);


% Define a ray as P + tu
% Where P is the point of origin of the ray, u is the unit vector direction
% of the ray, and t is the weight on that unit vector direction

ray = createLine3d(rayP0,theta,phi);
unitSphere = [0 0 0 1];

% Consider the transformation of the unit sphere to an ellipsoid by:
%  - Scaling in the p1, p2, and p3 directions
%  - Rotation about the p1, p2, and p3 directions by tip and tilt
%  - Translation
% We combine these in the 4x4 transformation matrix K.
% We then find the intersection of the ray inv(K)P + t*inv(K)u with the unit sphere.

% Construct K = translate * rotate * scale
scale = [ellipsoidRadii(1) 0 0 0; 0 ellipsoidRadii(2) 0 0; 0 0 ellipsoidRadii(3) 0; 0 0 0 1];
rotate = eye(4,4);
translate = eye(4,4);
translate(1:3,4) = ellipsoidCenter;
invK = (translate * rotate * scale)\eye(4,4);

% transform the ray
kP = invK*[ray(1:3) 0]';
ku = invK*[ray(4:6) 0]';
kRay = [kP(1:3)' ku(1:3)'];

% Find the intersection points
coordsOut = intersectLineSphere(kRay, unitSphere);

% Detect if we have a tangential ray
if coordsOut(1,:)==coordsOut(2,:)
    return
end

% Detect if we have a non-intersecting ray
if ~isnan(coordsOut(1,1))
    return
end

% If the radiiSign is positive, report the coordinates on the left-hand
% side of the ellipse, otherwise report the coordinates on the right
if radiiSign<0
    [~,rightIdx] = max(coordsOut(:,1));
    coordsOut = coordsOut(rightIdx,:);
else
    [~,leftIdx] = min(coordsOut(:,1));
    coordsOut = coordsOut(leftIdx,:);
end

% Now find the curvature of the ellipsoid at the intersection. To do so, we
% obtain the radii of the ellipse in the orthogonal planes of the rotated
% ray for the theta (p1p2) and phi (p1p3) angles.


% The radii of the ellipse that lies within the p1p2 plane when the ray is
% diverging from the optical axis into the p3 dimension by angle Phi
A = tan(phi); B = 0; C = 1;
D = coordsOut(2)/((coordsOut(1)-ellipsoidCenter(1))*A);
if isinf(D) || isnan(D)
    D = 0;
end
[ellipseRadii_p1p2(1),ellipseRadii_p1p2(2)]=EllipsoidPlaneIntersection(A,B,C,0,ellipsoidRadii(1),ellipsoidRadii(2),ellipsoidRadii(3));

% The radii of the ellipse that lies within the p1p3 plane when the ray is
% diverging from the optical axis into the p2 dimension by angle theta
A = tan(theta); B = 1; C = 0;
D = coordsOut(3)/((coordsOut(1)-ellipsoidCenter(1))*A);
if isinf(D) || isnan(D)
    D = 0;
end
[ellipseRadii_p1p3(1),ellipseRadii_p1p3(2)]=EllipsoidPlaneIntersection(A,B,C,0,ellipsoidRadii(1),ellipsoidRadii(2),ellipsoidRadii(3));

% Calculate the radius of curvature encountered by the ray.
% p1p2 plane
t = acos((coordsOut(1)-ellipsoidCenter(1))/ellipseRadii_p1p2(1));
curvature_p1p2 = radiiSign*((ellipseRadii_p1p2(1)^2*sin(t)^2 + ellipseRadii_p1p2(2)^2*cos(t)^2)^(3/2))/(ellipseRadii_p1p2(1)*ellipseRadii_p1p2(2));
% p1p3 plane
t = acos((coordsOut(1)-ellipsoidCenter(1))/ellipseRadii_p1p3(1));
curvature_p1p3 = radiiSign*((ellipseRadii_p1p3(1)^2*sin(t)^2 + ellipseRadii_p1p3(2)^2*cos(t)^2)^(3/2))/(ellipseRadii_p1p3(1)*ellipseRadii_p1p3(2));


end

