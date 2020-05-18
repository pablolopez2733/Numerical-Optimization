function [alpha] = zoom(alphalo, alphahi, f, xk, dk, gk)

%Definimos parametros y funciones utiles
c1 = 1e-4;
c2 = 0.99;
phid0 = dot(gk,dk);

phi = @(x) f(xk + x*dk);
L = @(y) f(xk) + c1*y*phid0;
phid = @(z) dot(apGrad(f,xk+z*dk), dk);

%Algoritmo zoom
while 1
    alphaj = (alphalo + alphahi)/2;
    
    if phi(alphaj) > L(alphaj)  || phi(alphaj) >= phi(alphalo)
        alphahi = alphaj;
    else
        if abs(phid(alphaj)) <= -c2*phid0
            alpha = alphaj;
            break
        end
    
        if phid(alphaj)*(alphahi-alphalo) >= 0
        alphahi = alphalo;
        end
        
        alphalo = alphaj;
        
    end
end

end

        
        