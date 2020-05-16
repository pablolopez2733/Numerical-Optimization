function [x, iter] = rcSR1( f, x0, itmax) 

% In : f ... (handle) function to be optimized 
%      x0 ... (vector) initial point 
%     itmax ... (natural number) upper bound for number of iterations 

% Out: x ... (vector) last approximation of a stationary point 
%     iter ... (natural number) number of iterations

%Step 1: We define key variables to be used
eta = 1e-2;
r = 1e-6;
tol = 1e-5;
rmax = 1.25;
radio = rmax;
iter = 0;
xk = x0;
n = length(xk);
g = apGrad(f,xk);
H = speye(n);
B = H;

%Usamos el pseudocódigo del Tema 4 Lab 2:

while norm(g, 'inf') > tol
    %P1) 
    s = -H*g;
    if dot(x,g) <0
        if norm(s) > radio
            s = radio*s/norm(s);
        end
    else 
        s = pCaucht(B,g,radio)
    end
    %P2) Calcular la reducción del modelo y de la función:
    new = xk + s;
    redf = f(xk) - f(new);
    redm = - (dot(g,s)+0.5*dot(s,B*s));
    div = redf/redm;
    gnew = apGrad(f, new);
    gamma = gnew-g;
    
    %P3) Si redf/redm > η, entonces xk+1 def == xk + sk . En otro caso xk+1 def == xk.
    if (div > eta)
        xk = xk + s;
        g = gnew;
        iter = iter + 1;
        
        if iter == itermax
            break
        end
    end
    
    % P4) Si redf/redm > 0.75, entonces Si ||sk|| > 0.8∆k, entonces ∆k+1 def== 2∆k. En otro caso ∆k+1 def== ∆k.
    if div > 0.75
        if norm(s)>0.8*radio;
            radio = min(2*radio,rmax);
        end
    %P5) Si redf/redm < 0.1, entonces ∆k+1 def == 0.5∆k    
    elseif div < 0.1
        radio = 0.5*radio;
    end
    %P6)  Si se cumple la condici´on (1), entonces actualizamos las matrices Bk+1 y Hk+1 usando las formulas (SR1) y (iSR1). Luego asignar k ← k + 1 y continuar con el ciclo. 
    v = gamma - B*s;
    if abs(dot(v,s))>=r*norm(s)*norm(v)
        B = (v/dot(v,s))*v' + B;
        u = s - H*gamma;
        H = (u/dot(u,gamma))*u' + H;
    end
end

        
    


