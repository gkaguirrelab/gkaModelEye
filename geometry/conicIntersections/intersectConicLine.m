% intersectConicLine - given a line and a conic detect the real intersection
% points 
%          - Pierluigi Taddei (pierluigi.taddei@polimi.it)
%
% Usage:   P = intersectConicLine(C, l)
%
% Arguments:
%           C - homogeneous conic matrix
%           l - homogeneous line vector
%
%           P - matrix of homogeneous intersection points (each column is a
%           point)
% 08.3.2007 : Created
%
function P = intersectConicLine(C, l)
    P = [];
    [p1 p2] = getPointsOnLine(l);  
    
    %p1'Cp1 + 2k p1'Cp2 + k^2 p2'Cp2 = 0
    p1Cp1 = p1' *C* p1;
    p2Cp2 = p2' *C* p2;
    p1Cp2 = p1' *C* p2;

    if (p2Cp2 == 0) %linear
       k1 = -0.5*p1Cp1 / p1Cp2;
       P = [p1 + k1*p2];
    else
        delta = p1Cp2^2 - p1Cp1*p2Cp2;
        if (delta >= 0)
            deltaSqrt = sqrt(delta);
            k1 = (-p1Cp2 + deltaSqrt)/p2Cp2;
            k2 = (-p1Cp2 - deltaSqrt)/p2Cp2;

            P = [p1 + k1*p2, p1 + k2*p2];
       % else
       %     disp('no intersection points detected');
        end

    end
%     %the intersection is managed as the decomposition of the dual conic
%     Ml = crossMatrix(l);
%     B = Ml'*(l'*l)*Ml;
%     
%     if (l(3) ~= 0)
%         a = sqrt(C(1,1)*C(2,2)- C(1,2)^2)/l(3);
%     else
%         disp('infinite line are not covered: implements this branch...');
%     end
%     
%     D = B + a*Ml;
%     
%     ci = find(D ~= 0);
%     j = floor(ci(1) / 3)+1;
%     i = ci(1) - (j-1)*3;
% 
%     p = D(i,:)';
%     q = D(:,j);
%     
%     P = [];
%     if (isreal(p)) 
%         P = p;
%     end
%    	if (isreal(q)) 
%         P = [P, q];
%     end

end