function [x, msg, iter] = mRCH(f, x0, itmax)
% Trust region method using the dogleg point
% 
% In:   f     ... (handle) function to be optimized
%       x0    ... (vector) initial point
%       itmax ... (natural number) upper bound for number of iterations
% 
% Out:  x     ... (vector) last approximation of a stationary point
%       msg   ... (string) message that says whether (or not) a minimum 
%                          was found
    xk = x0;
    delta = 0.5;
    eta = 0.1;
    tol = 1e-5;
    deltaMax = 1.5;
    n = length(x0);
    
    for k = 1:itmax
        fk = f(xk);
        gk = apGrad(f, xk);
        Bk = apHess(f, xk);
        
        lambda1 = eigs(Bk, 1, 'smallestreal');
        s = 1e-12 * - 9/8 * lambda1;
        
        if lambda1 <= 0
            Bk = Bk + s * eye(n);
        end 
        
        mk = @(p) fk + dot(gk, p) + 0.5 * dot(p, Bk*p);
        pk = pDogLegH(Bk, gk, delta);
        
        rhok = (f(xk) - f(xk + pk)) / (mk(zeros(n, 1)) - mk(pk));
        
        if rhok < 0.25
            delta = 0.25 * delta;
        elseif rhok > 0.75 && norm(pk) == delta
            delta = min(2 * delta, deltaMax);
        end

        if rhok > eta
            xk = xk + pk;
        end

        if norm(gk) < tol
            x = xk;
            msg = 'Se encontro el minimo :)';
            iter = k;
            break
        elseif k == itmax
            x = xk;
            iter = itmax;
            msg = 'No se encontro el minimo :(';
        end
    end
    
end