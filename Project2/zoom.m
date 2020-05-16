function [alphaStar] = zoom(alphaLow, alphaHigh, f, xk, gk, dk)

%Define Parameters
c1 = 1e-4;
c2 = 0.99;

%Define Functions 
slope = dot(gk,dk);
phi = @(x) f(xk + x*dk);
L = @(y) f(xk) + c1*y*slope;
phiP = @ (z) dot(apGrad(f,xk+z*dk), dk);

while 1
    alphaMid = (alphaLow + alphaHigh)/2;
    if phi(alphaMid)>=phi(alphaLow) || phi(alphaMid)> L(alphaMid)
        alphaHigh = alphaMid;
    elseif abs(phiP(alphaMid))<=-c2*slope
        alpha2 = alphaMid;
        break
    elseif phiP(alphaMid)*(alphaHigh-alphaLow)>=0
        alphaHigh = alphaLow;
        alphaLow = alphaMid;
    else 
        alphaLow = alphaMid;
    end
end
alphaStar = alpha2;
end

        
        