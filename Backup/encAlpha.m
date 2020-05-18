function [alpha, gnew] = encAlpha( f, xk, dk, gk )
% Purpose: find the first alpha of the form (1/2)^k that satisfies the
% If (W2) is not satisfied, throw an error.
%
% In : f ... function to minimize
%      xk ... current point
%      dk ... descend direction
%  gk ... gradient of f in xk
%
% Out: alpha ... parameter satisfying the two Wolfe conditions.
%gnew ... gradient of f in xk + alpha*dk
%
% parameters: c1 = 1e-4, c2 = 0.99,alpha0 = 1
%

    c1 = 1e-4;
    c2 = 0.99;
    alpha = 1;
    
    fk = f(xk);
    slope = dot(dk, gk);
    while f(xk + alpha*dk) > fk + c1*alpha*slope
        alpha = 0.5*alpha;
    end
    gnew = apGrad(f, xk + alpha*dk);
    assert( dot(gnew, dk) >= c2*slope, '(W2) no se cumple')
end