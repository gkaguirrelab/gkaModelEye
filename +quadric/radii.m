function r = radii(S)


% If the quadric surface was passed in vector form, convert to matrix
if isequal(size(S),[1 10])
    S = quadric.vecToMatrix(S);
end

% find the center of the quadric
center = -S( 1:3, 1:3 ) \ S( 1:3,4 );

% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';

% translate to the center
Q = T * S * transpose(T);

% solve the eigenproblem
[~,evals] = svd(Q( 1:3, 1:3 ) / -Q( 4, 4 ));
r = 1./diag(sqrt(evals));
sgns = sign( diag( evals ) );
r = r .* sgns;

% We now have the radii. We would like to return these values in the order
% of the dimensions of the space (x, y, z). As a convention, the radius
% along the (e.g.) x dimension is the radius that correspond to the axis of
% the quadric surface that has an angle between -45 and 45 degrees w.r.t.
% to the x axis.

% Get the orientation
thisOrientation = orientFxn(S);

% % Define a mapping of orientations to radii order
orientations = {[0 0 0],[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[1 1 1]};
orders = {[1 2 3],[1 3 2],[3 2 1],[2 1 3],[2 3 1],[3 1 2],[1 2 3]};
thisOrder = orders{cellfun(@(x) isequal(x,thisOrientation),orientations)};

% Report the radii in the specified order
r = r(thisOrder);

end


function orientation = orientFxn(S)
a = quadric.angles(S);
orientation = mod(fix((abs(a)+45)./90),2);
end