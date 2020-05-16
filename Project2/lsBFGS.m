function [x, iter] = lsBFGS( f, x0, itmax)
% In : f ... (handle) function to be optimized 
%      x0 ... (vector) initial point 
%     itmax ... (natural number) upper bound for number of iterations 

% Out: x ... (vector) last approximation of a stationary point 
%     iter ... (natural number) number of iterations

%1. We define key variables
iter = itmax;
xk = x0;
n = length(x0);
H = speye(n);
g = apGrad(f,xk);
tol = 1e-5

for k = 1:itmax
    if norm(g,'inf')<= tol
        iter = k-1;
        break
    end
    dk = -H*g;
    [alpha,gnew] = lineSearch(f,xk,dk,g);
    s = alpha*dk;
    xk = xk + s;
    gamma = gnew - g;
    rhoinv = dot(gamma,s);
    Hgamma = H*gamma/rhoinv;
    Hgammas = s*Hgamma' + Hgamma*s';
    H = H - Hgams + ((1+dot(gamma,Hgamma)/(rhoinv*s))*s';
    g = gnew;
end
    xf = xk;
end


   

