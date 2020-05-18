function [xk, iter] = rcSCR1( f, x0, itmax) 

% In : f    ... (handle) function to be optimized 
%      x0   ... (vector) initial point 
%     itmax ... (natural number) upper bound for number of iterations 

% Out: x ... (vector) last approximation of a stationary point 
%     iter ... (natural number) number of iterations

%Definimos parametros y funciones utiles
eta = 0.1;
r = 1e-6;
tol = 1e-5;
deltamax = 1.25;
deltak=deltamax;
n = length(x0);
xk=x0;
Hk = speye(n);
Bk = Hk;
iter = 0;
gk=apGrad(f,x0);

%Metodo SR1 Trust-Region 
while norm(gk) > tol && iter < itmax
    sk = -Hk*gk;
    if dot(sk,gk) < 0
        if norm(sk) > deltak
            sk = deltak*sk/norm(sk);
        end
    else 
        sk = pCauchy(Bk,gk,deltak);
    end
    
    yk = apGrad(f,xk+sk) - gk;
    redf = f(xk) - f(xk+sk);
    redm = - (dot(gk,sk)+0.5*dot(sk,Bk*sk));
    div = redf/redm;
    
    if div > eta
        xk = xk + sk;
        gk = apGrad(f,xk);
        iter = iter + 1;    
    end
    
    if div > 0.75
        if norm(sk)>0.8*deltak
            deltak = 2*deltak;
        end
    end    
    if div < 0.1
        deltak = 0.5*deltak;
    end
    
    vk = yk - Bk*sk;
    
    if abs(dot(vk,sk))>=r*norm(sk)*norm(vk)
        Bk = Bk + (vk*vk')/dot(vk,sk);
        wk=sk-Hk*yk;
        Hk = Hk + (wk*wk')/dot(wk,yk);
    end
end

%if norm(gk)<= tol
    %x=xk;
%end
end


        
    


