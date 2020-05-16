function [alphastar, gnew] = lineSearch( f, xk, dk, gk ) 
% In : f ... objectve function (handle) 
%      xk ... current point 
%      dk ... chosen direction of descent 
%      gk ... gradient of f in xk % 
% Out: alpha ... a parameter satisfying (W1) and (W2) 
%      gnew ... gradient of f in xk+alpha*dk

%Define Parameters
alphaold = 0;
alphanew = 1;
alphamax = 2^10;
c1 = 1e-4;
c2 = 0.99;

%Define Useful Functions
slope = dot(gk,dk);
phi = @(x) f(xk + x*dk);
L = @(y) f(xk) + c1*y*slope;
phiP = @ (z) dot(apGrad(f,xk+z*dk), dk);

%LineSearch
while alphanew > 0 && alphanew<alphamax
    if phi(alphanew)>leg(alphanew) || phi(alphanew) >= phi(alphaold)
        alpha2 = zoom(alphaold,alphanew,f,xk,gk,dk);
        break
    elseif abs(phiP(alphanew))<=-c2*slope
        alpha2 = alphanew;
    break
    elseif phiP(alphanew)>=0
        alpha2 = zoom(alphanew, alphaold,f,xk,gk,dk);
    break
    else
        alphaold = alphanew;
        alphanew = min(2*alphanew, alphamax - 10*eps);
    end
end
alphastar = alpha2;
gnew = apGrad(f,xk+alphastar*dk);
