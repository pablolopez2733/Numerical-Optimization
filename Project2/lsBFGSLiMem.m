function [x, iter] = lsBFGSLiMem( f, x0, itmax, m) 
% In : f ... (handle) function to be optimized 
%      x0 ... (vector) initial point 
%      itmax ... (natural number) upper bound for number of iterations 
%      m ... (natural number) number of recent update steps (limited memory parameter)
% Out: x ... (vector) last approximation of a stationary point 
%     iter ... (natural number) number of iterations

%1. We define key variables
iter = itmax;
xk = x0;
n = length(x0);
H = speye(n);
g = apGrad(f,xk);
S = zeros(n,m);
tol = 1e-5
Gamma = S;
dk = -g;

for k = 1:itmax
    if norm(g,'inf')<= tol
        iter = k-1;
        break
    end
    assert(dot(dk,g)<0)
    [alpha, gnew] = lineSearch(f,xk,dk,g);
    
    s = alpha*dk;
    xk = xk + s;
    S = [alpha*dk,S(:,1:m-1)];
    Gam = [gnew - g, Gam(:,1:m-1)];
    
    if k<=m
        dk = -findHg(S(:,1:k),Gamma(:,1:k),gnew);
    else
        dk = -findHg(S,Gamma,gnew);
    end
    g = gnew;
end
end


