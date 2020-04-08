


function [p] = pDogLeg(B,g,delta)
%Our code for this is based on the flowchart that we received in class.
%We compute alpha u:
alphau=(g.'*g)/(g.'*B*g);
%Check if it is greater than delta over norm of g
if alphau >= delta/norm(g)
    p=(-1)*(delta/norm(g))*g;%The point is the Cauchy Point
else
    Binv=inv(B);
    pb=(-1)*Binv*g;
    
    if norm(pb)<= delta
        p=pb;%The point is Newton´s point
    else
        %We compute Pu:
        pu= (-.99)*alphau*g;
        %Now we have to find alpha star between 0 and 1
        %by solving: ||Pu+alpha(Pb-Pu)||^2
        coef_alpha=[(pb-pu).'*(pb-pu),2*(pu.'*pb-pu.'*pu),(pu.' * pu)-delta^2];
        sols=roots(coef_alpha);
        if(sols(1)>0 && sols(1)<=1)
            alpha_star=sols(1);
        else
            alpha_star=sols(2);
        end
        
        
        p=pu+(alpha_star)*(pb-pu);
        
        
    end
end
        

end