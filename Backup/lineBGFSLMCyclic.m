function [xk, iter] = lineBGFSLMCyclic( f, x0, tol, maxiter, m)

%Este método es una modificación al método lineBGFSLM realizado en clase
%que incorpora memoria cíclica.  
%La modificación consiste en que ahora, en cada iteración, guardamos los
%nuevos elementos s y gnew en una posición variable en lugar de al
%principio de las matrices S y Gam (respectivamente). La posición en la que
%se guardan los nuevos elementos está dada por la variable pointer, que
%aumenta en cada iteración y, cuando llega a m, regresa a la posición
%inicial pointer = 1. 

% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
%
% In : f ... function to minimize
% x0 ... initial point
% tol ... tolarance
% maxiter ... upper bound for iterations
% m ... memoria (how much history)
%
% Out: xf ... final approximation of x*
% iter ... number of iterations used
%
% Initial values: H = I
% Criterios de paro: || gk || <= tol or ||s|| <= 1e-7
%

    n = length(x0);
    iter = maxiter;
    xk = x0;
    g = apGrad( f, xk );
    %H = speye(n);
    S = zeros(n, m);   Gam = S;

    dk = -g;
    pointer = 1;
   for k = 1:maxiter
       if norm(g, 'inf') <= tol
           iter = k-1;
           break
       end
       assert( dot(dk, g) < 0)
       
       [alpha, gnew] = lineSearch(f, xk, dk, g ); %P2
       
       %memorizar
       s   = alpha*dk;
       xk  = xk + s;
       S(:,pointer) = s;
       Gam(:,pointer) = gnew-g; 
       
       if k <= m
           dk = -calcHg( S(:, 1:pointer), Gam(:, 1:pointer), gnew);
       else
           dk = -calcHg([S(:,pointer:m),S(:,1:pointer-1)], [Gam(:,pointer:m),Gam(:,1:pointer-1)], gnew);
       end
       if (m ~= 1)
           if(pointer == m-1)
               pointer = m;
           else
               pointer = mod(pointer + 1,m);
           end
       end
       g = gnew;
   end
end


function [q] = calcHg(S, Gam, gnew )
    % S = {mas nuevo ... mas viejo}
    q = gnew;
    m = length(S(1,:));
    irhos = dot(S, Gam, 1);
    alphas = zeros(m,1);
    
    for i = 1:m
        alphas(i) = dot(S(:,i), q)/irhos(i);
        q = q - alphas(i)*Gam(:,i);
    end
    
    delta = irhos(1)/dot(Gam(:,1), Gam(:,1));
    q = delta*q;
    
    for i = m:-1:1
       beta = dot(Gam(:,i), q)/irhos(i);
       q = q +(alphas(i) -beta)*S(:,i);
    end 
end