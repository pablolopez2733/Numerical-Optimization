function [x,msg] = mRC1(f,x0,itmax)
% Trust region method using the dogleg point.
% Same arguments and results as the mRC1 method:
% In : f ... (handle) function to be optimized
% x0 ... (vector) initial point
% itmax ... (natural number) upper bound for number of iterations
%
% Out: x ... (vector) last approximation of a stationary point
% msg ... (string) message that says whether (or not) a minimum was found

%For both RC1 and RC2 we based our code on the Algorithm 4.1 of Nocedal�s
%Numerical Optimization book.

del0=1; %Delta 0
delm=1.5; %Delta Max
eta=0.10; %eta
tol=0.00001;
format long E;

g0=apGrad(f, x0);
mg0=norm(g0);

if (mg0<=tol)
  x=x0;
  msg="El m�nimo fue hallado";

else
  mgk=mg0;
  delk=del0;
  xk=x0;
  k=1;

  while (k<itmax && mgk>tol)

    gk=apGrad(f,xk);
    mgk=norm(gk);
    Bk=apHess(f,xk);
    pck=pCauchy(Bk,gk,delk);
    
    rhok=-(f(xk)-f(xk+pck))/(dot(gk,pck)+0.5*(dot(pck,Bk*pck)));
        
    mindel=min(2*delk,delm);
        
    if (rhok<0.25) 
      delk1=0.25*delk;
    else
      if ((rhok>0.75)&&(norm(pck)==delk))
        delk1=mindel;
      else
        delk1=delk;
      end
    end
    
    if (rhok>eta)
       xk1=xk+pck;
    else
       xk1=xk;     
    end
    
    gk1=apGrad(f,xk1);
    mgk1=norm(gk1);
    mgk=mgk1;
    delk=delk1;
    xk=xk1; 
    k=k+1;
  
  end
x=xk;
if mgk<tol
    msg="El m�nimo fue hallado";
else
    msg="El m�nimo NO fue hallado";
end


end   

end
