function [xk, iter] = lsBFGS( f, x0, itmax)
% In : f ... (handle) function to be optimized 
%      x0 ... (vector) initial point 
%     itmax ... (natural number) upper bound for number of iterations 

% Out: xk ... (vector) last approximation of a stationary point 
%     iter ... (natural number) number of iterations

%1. We define key variables
iter = 0;
xk = x0;
n = length(x0);
H = speye(n);
gk = apGrad(f,xk);
tol = 1e-5;
gnew = 1;

    while norm(gnew) > tol
        dk = -H*gk;
    
        [alpha, gnew] = encAlpha( f, xk, dk, gk );
    
        s = alpha*dk;
        xk = xk + s;
        gamma = gnew - gk;
    
        rhoinv = dot(s, gamma);
        assert(rhoinv > 0);
    
        Hgamma = H*gamma/rhoinv;
        H = H -(s*Hgamma' + Hgamma*s') + ((dot(gamma, Hgamma)+1)/rhoinv*s)*s';
        gk = gnew;
        iter = iter+1;   % STOP Condition
        if  iter >= itmax
            break
        end
    end
end