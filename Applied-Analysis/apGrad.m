function [g] = apGrad(f, x)
% In : f ... (handle) anonymous function
%      x ... (vector) point where to approx. the gradient

    n = length(x);
    g = zeros(n,1);
    %hs = nthroot(eps,3)*(speye(n) + spdiags(abs(x), 0, n,n));
    hs = nthroot(eps,3)*spdiags(1+abs(x), 0, n,n);
    
    for i = 1:n
        g(i) = 0.5*( f(x + hs(:,i))-f(x - hs(:,i)) )/hs(i,i);
    end
end