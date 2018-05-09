% crossMatrix - given a homogeneous point generates the skew symmetric
% matrix form used in cross products: p x q -> Mp . q
%          - Pierluigi Taddei (pierluigi.taddei@polimi.it)
%
% Usage:   Mp = adjoint3(M)
%
% Arguments:
%           p - homogeneous 2d point
%           Mp - 3x3 skew symmetric matrix
% 08.3.2007 : Created
%
function Mp = crossMatrix(p)
    Mp = zeros(3,3);
    Mp(1,2) = p(3);
    Mp(1,3) = -p(2);
    Mp(2,1) = -p(3);
    Mp(2,3) = p(1);
    Mp(3,1) = p(2);
    Mp(3,2) = -p(1);
end