function [H] = apHess(f, x)
% In : f ... (handle) anonymous function
%      x ... (vector) point where to approx. the gradient

    n = length(x);
    H = zeros(n,n);
    hs = nthroot(eps,4)*spdiags(1+abs(x), 0, n,n);
    
    for j = 1:n
    for i = j+1:n
        H(i, j) = 0.25*( f(x + hs(:,i) + hs(:,j)) ...
                       + f(x - hs(:,i) - hs(:,j)) ...
                       - f(x + hs(:,i) - hs(:,j)) ...
                       - f(x - hs(:,i) + hs(    :,j)) )...
                       /( hs(i,i)*hs(j,j) );
    end
    end
    
    H = H + H';
    
    for j = 1:n
        H(j,j) = ( f(x + hs(:,j)) + f(x - hs(:,j)) - 2*f(x) )...
                       /( hs(j,j)*hs(j,j) );
    end    
end