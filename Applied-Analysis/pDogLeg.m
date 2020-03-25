


function [p] = pDogLeg(B,g,delta)
%Vamos a usar el diagrama de flujo visto en clase
%Calculamos alpha u:
alphau=(g.'*g)/(g.'*B*g);
%Checamos si es mayor a delta entre norma de g:
if alphau >= delta/norm(g)
    p=(-1)*(delta/norm(g))*g;%Punto es pto Cauchy
else
    pb=(-1)*inv(B)*g;
    
    if norm(pb)<= delta
        p=pb;%punto es punto de Newton
    else
        %Calculamos entonces Pu:
        pu= (-.99)*alphau*g;
        %Ahora hay que encontrar alpha asterisco entre 0 y 1
        %Resolviendo: ||Pu+alpha(Pb-Pu)||^2
        coef_alpha=[(pb-pu).'*(pb-pu),2*(pu.'*pb-pu.'*pu),(pu.' * pu)-delta^2];
        alpha_star=roots(coef_alpha);
        p=pu+(alpha_star)*(pb-pu);
        
        
    end
end
        

end