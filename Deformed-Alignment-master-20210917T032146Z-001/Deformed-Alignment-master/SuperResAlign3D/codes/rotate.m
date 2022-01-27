 function [xnew,ynew]=rotate(theta0,xold,yold)

    xnew = xold * cos(theta0) - yold * sin(theta0);
    ynew = xold * sin(theta0) + yold * cos(theta0);
    
 % plot(xnew,ynew,'.b'); axis equal;
 hold on