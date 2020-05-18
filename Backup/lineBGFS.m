function [xk, k] = lineBGFS( f, x0, tol, maxiter )
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
%
% Same parameters and results as lineDFP,
% with the exception that the (iBGFS) update formula is used.
%

    k = 0;
    n = length(x0);
    I = speye(n);
    H = speye(n);
    gk = apGrad(f, x0);
    gnew = 1;
    xk = x0;
    while norm(gnew) > tol
        dk = -H*gk;
        
        [alpha, gnew] = encAlpha( f, xk, dk, gk );
        
        s = alpha*dk;
        xk = xk + s;
        gamma = gnew - gk;
        rhoinv = dot(s, gamma);
        
        H = (I - s*gamma'/rhoinv)*H*(I - gamma*s'/rhoinv) + (s*s')/rhoinv;  % O(n^3)
        gk = gnew;
        k = k+1;   % STOP si es grade.
        if  k >= maxiter
            break
        end
    end
end