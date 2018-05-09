% PLOTCONIC: Plots the matrixr form of an arbitrary conic:
%            Ax^2 + By^2 + 2Cx + 2Dy + Exy + F = 0
%
%        Conservative upper and lower bounds on x can be specified,
%        and will be adjusted inward.
%
%     Usage: conic(K,xBounds)
%
%             K = conic matrix
%             xBounds, yBounds = maximum x upper and lower bounds, by default [-50,50]
%
% the homogeneous representation of a conic is a matrix
% m = [A C D; C B E; D E F] that represents the equation
% A x^2 + B y^2 + 2C xy + 2D x + 2Ey + F = 0


function plotConic(K,color,xBounds, yBounds)

    K = K / norm(K); %normalize
    a = K(1,1);
    b = K(2,2);
    e = K(1,2)*2;
    c = K(1,3)*2;
    d = K(2,3)*2;
    f = K(3,3);

    %the equation is ax^2 + by^2 + cx + dy + exy + f = 0
    if (nargin < 2)
       color = 'b';
       xBounds = [-50,50];
       yBounds = [-50,50];
    elseif (nargin < 3)
        xBounds = [-50,50];
        yBounds = [-50,50];
   elseif (nargin < 4)
        yBounds = xBounds;
    end

  x = linspace(xBounds(1),xBounds(2),5000)';
  [ypos yneg] = solveQuadratic(a,b,c,d,e,f,x);

  i = find(( ~imag(ypos) | ~imag(yneg) ));    % Find bounds of soln
  leni = length(i);
  if (leni == 0)
    disp('  No real solution found');
  else
  
    xBounds(1) = x(max([min(i),1]));
    xBounds(2) = x(min([max(i),length(x)]));
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SOlVE SYSTEMS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

    
    x = linspace(xBounds(1), xBounds(2),200)';  % New x bounds
    [ypos yneg] = solveQuadratic(a,b,c,d,e,f,x);
  
    iypos = ~imag(ypos) & ypos > yBounds(1) & ypos < yBounds(2);
    iyneg = ~imag(yneg) & yneg > yBounds(1) & yneg < yBounds(2);

    good = [yneg(iyneg); ypos(iypos)];
    if (isempty(good))
         disp('No solution for current bounds');
         return;
    end
    
    plot(x(iypos),ypos(iypos), color);       % Plot pos y as fn(x)
    hold on;
    plot(x(iyneg),yneg(iyneg), color);       % Plot neg y as fn(x)
    
    upper = max(good);
    lower = min(good);

    y = linspace(lower, upper, 200)'; % y bounds

    
    [xpos xneg] = solveQuadratic(b,a,d,c,e,f,y);
    ixpos = ~imag(xpos) & xpos > xBounds(1) & xpos < xBounds(2);
    ixneg = ~imag(xneg) & xneg > xBounds(1) & xneg < xBounds(2);



    plot(xpos(ixpos),y(ixpos), color);       % Plot pos x as fn(y)
    %config = [config; xpos(ixpos) y(ixpos)];

    plot(xneg(ixneg),y(ixneg), color);       % Plot neg x as fn(y)
    %config = [config; xneg(ixneg) y(ixneg)];

    %bounds = sqplot(config(:,1),config(:,2));
    %axis(bounds);
    %axis('equal');
    %hold off;
  end;

  return;
end

function [yPos yNeg] = solveQuadratic(a,b,c,d,e,f, x)
  ka = b*ones(size(x));      % Quadratic coefficients of y
  kb = (e.*x + d) ;      % First degree coefficents of y
  kc = a*(x.^2) + c*x + f;    
  
if (b == 0) %%linear system
    
    yPos = -kc./kb;
    yNeg = yPos;
    
else

  % Solve quadratic
  delta = sqrt( kb.^2 - (4.* ka.* kc));
  yPos = (-kb + delta) ./ (2*ka);
  yNeg = (-kb - delta) ./ (2*ka);
  
end
  %test!
  x = [x;x];
  y = [yNeg; yPos];
  
  te = a * x.^2 + b * y.^2 + c.*x + d.*y + e.*x.*y + f;
  max(te)
  
end