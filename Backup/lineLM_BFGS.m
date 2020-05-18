function [xk, k] = lineLM_BFGS( f, x0, tol, itmax,m )
% Purpose: use limited memory updates to iapproximate the inverse of the 
%          Hessian. (linesearch algorithm)
%

k = 0;
SS = [];
GG = [];
gk = apGrad(f, x0);
gnew = 1;
xk = x0;
    while norm(gnew) > tol
        dk = -evaluate_Hg(SS, GG, gk);
        
        [alpha, gnew] = lineSearch( f, xk, dk, gk );
        
        s = alpha*dk;
        xk = xk + s;
        gamma = gnew - gk;
        %rhoinv = dot(s, gamma);
        
        % we have to refresh the memory
        if k <= m
            SS = [s, SS];
            GG = [gamma, GG];
        else
            SS = [s, SS(:, 1:m-1)];
            GG = [gamma, GG(:, 1:m-1)];
        end
        
        gk = gnew;
        k = k+1;   % STOP si es grade.
        if  k >= itmax
            break
        end
    end
end


function r = evaluate_Hg(SS, GG, gk)
    % SS ... s(k-1), s(k-2), ... s(k-m)
    % GG ... gam(k-1), gam(k-2), ... gam(k-m)
    % gk ... gradiente en xk.
    
    [~, m] = size(GG);
    
    if m == 0
        r = gk;
    else
        alphas = zeros(m,1);
        r = gk;
        for j = 1:m
            alphas(j) = dot( SS(:,j), r )/dot( SS(:,j),GG(:,j) );
            r = r - alphas(j)*GG(:,j);
        end
        
        r = dot( SS(:,1),GG(:,1) )/dot( GG(:,1),GG(:,1) )*r;
        
        for j = m:-1:1
            beta = dot( GG(:,j), r )/dot( SS(:,j),GG(:,j) );
            r = r + (alphas(j)-beta)*SS(:,j);
        end
    end
end