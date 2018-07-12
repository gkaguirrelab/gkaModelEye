function axisOrder = axisOrder(S)

% We now have the radii. We would like to return these values in the order
% of the dimensions of the space (x, y, z). As a convention, the radius
% along the (e.g.) x dimension is the radius that correspond to the axis of
% the quadric surface that has an angle between -45 and 45 degrees w.r.t.
% to the x axis.
%{
    S = quadric.unitSphere;
    S = quadric.scale(S,[3 5 4]);
  r = quadric.radii(S)
  axisOrder = axisOrientation(S);
% Report the radii in the specified order
r = r(axisOrder)
%}

% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

a = quadric.angles(S);
thisOrientation = mod(fix((abs(a)+45)./90),2);

% Define a mapping of orientations to radii order
orientations = {[0 0 0],[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[1 1 1]};
orders = {[1 2 3],[1 3 2],[3 2 1],[2 1 3],[2 3 1],[3 1 2],[1 2 3]};
axisOrder = orders{cellfun(@(x) isequal(x,thisOrientation),orientations)};


end