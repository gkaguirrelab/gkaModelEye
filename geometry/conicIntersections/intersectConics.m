% intersectConics - intersects two non degenerate conics
%          - Pierluigi Taddei (pierluigi.taddei@polimi.it)
%
% Usage:   P = intersectConics(E1, E2)
%
% Arguments:
%           E1, E2 - homogeneous conic matrices. E1 must be non degenerate
%
%           P - matrix of homogeneous intersection points (each column is a
%           point)
% 08.3.2007 : Created
%
function P = intersectConics(E1, E2)
 
    r1 = rank(E1);
    r2 = rank(E2);
    
    if (r1 == 3 && r2 == 3)
        P = completeIntersection(E1,E2);
        return
    else
        if (r2 < 3) %E2 is degenerate
            defE = E2;
            fullE = E1;
        else
            defE = E1; %E1 is degenerate
            fullE = E2;
        end
        
        [m l] = decomposeDegenerateConic(defE);
        P1 = intersectConicLine(fullE,m);
        P2 = intersectConicLine(fullE,l);
        P = [P1, P2];
    end
end

function P = completeIntersection(E1,E2)
 %exploit determinant multilinearity 
 %   k = [   det(E1); ...
 %           det([E1(:,1), E1(:,2), E2(:,3)]) + det([E1(:,1), E2(:,2), E1(:,3)]) + det([E2(:,1), E1(:,2), E1(:,3)]); ...
 %           det([E1(:,1), E2(:,2), E2(:,3)]) + det([E2(:,1), E1(:,2), E2(:,3)]) + det([E2(:,1), E2(:,2), E1(:,3)]); ...
 %           det(E2); ...
 %   ];

 %characteristic polynom: C1*(-C2)^-1 - lambda I
    EE = E1 * inv(-E2);
    k = [  -1; ...
             trace(EE);
           - ( det(EE(1:2,1:2)) + det(EE(2:3,2:3)) + det(EE([1,3],[1,3])) ); ...
            det(EE); ...
    ];
    
    r = roots(k);
    m = [];
    if (isreal(r(1))) 
        E0 = E1 + r(1)*E2;
        [m l] = decomposeDegenerateConic(E0);
    end
    if (isempty(m) && isreal(r(2)))
        E0 = E1 + r(2)*E2;
        [m l] = decomposeDegenerateConic(E0);
    end   
    if (isempty(m) && isreal(r(3)))
        E0 = E1 + r(3)*E2;
        [m l] = decomposeDegenerateConic(E0);
    end

     if (isempty(m))
      %  disp('no intersecting lines detected');
        P = [];
        return;
     end
    
    P1 = intersectConicLine(E1,m);
    P2 = intersectConicLine(E1,l);
    P = [P1, P2];
end