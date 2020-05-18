function [alpha, gnew] = lineSearch( f, xk, dk, gk ) 
% In : f  ... (handle) objectve function 
%      xk ... current point 
%      dk ... chosen direction of descent 
%      gk ... gradient of f in xk
%
% Out: alpha ... a parameter satisfying (W1) and (W2) 
%      gnew  ... gradient of f in xk+alpha*dk

%Definimos parametros y funciones útiles
alpha0 = 0;
alpha1 = 1;
alphamax = 1e3;
c1 = 1e-4;
c2 = 0.99;
phid0 = dot(gk,dk);
n=length(xk);

phi = @(x) f(xk + x*dk);
L = @(y) f(xk) + c1*y*phid0;
phid = @(z) dot(apGrad(f,xk+z*dk), dk);

%Algoritmo LineSearch
while alpha1>0 && alpha1<alphamax
    
    if phi(alpha1) > L(alpha1) || phi(alpha1) >= phi(alpha0)
        alpha = zoom(alpha0,alpha1,f,xk,gk,dk);
        gnew = apGrad(f,xk+alpha*dk);
        break
    end
    
    if abs(phid(alpha1)) <= -c2*phid0
        alpha = alpha1;
        gnew = apGrad(f,xk+alpha*dk);
        break
    end
    
    if phid(alpha1) >= 0
        alpha = zoom(alpha1, alpha0,f,xk,gk,dk);
        gnew = apGrad(f,xk+alpha*dk);
        break
    end
    
    alpha0 = alpha1;
    alpha1 = 2*alpha1;   
end

if alpha1 >= alphamax
    alpha=alpha1;
    gnew=apGrad(f,xk+alpha*dk);
end

end
