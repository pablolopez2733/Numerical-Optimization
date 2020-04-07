function [x,msg] = mRC2(f,x0,itmax)
% Trust region method using the dogleg point.
% Same arguments and results as the mRC1 method:
% In : f ... (handle) function to be optimized
% x0 ... (vector) initial point
% itmax ... (natural number) upper bound for number of iterations
%
% Out: x ... (vector) last approximation of a stationary point
% msg ... (string) message that says whether (or not) a minimum was found

%Define central parameters:
eta=.1;
tol=1e-5;
deltaMax=1.5;
delta=0+1.5*rand(1,1);
xk=x0;
k=1;
bandera=false;
while (k <= itmax)
    %Define values for each iteration:
    fk=f(xk);
    g=apGrad(f,xk);
    B=apHess(f,xk);
    %Check if B is spd
    lamda1=eigs(B,1,'sm');
    %If it�s not, we redefine B as follows:
    if(lamda1<=0)
        s=(1e-12)-((9/8)*(lamda1));
        B=B+s*(eye(length(xk)));
        
    end
    %Following the Algorithm 4.1 on Nocedal,we first Obtain pk by
    %approximately solving: mk(p) s.t. ||p||<=delta.
    %This is where we use our dogleg algorithm to approximate p
    mk= @(p) fk+g.'*p+(1/2)*p.'*B*p;
    p=pDogLeg(B,g,delta);
    %Evaluate rho from 4.4
    rho=(f(xk)-f(xk+p))/(mk(zeros(size(xk)))-mk(p));
    if (rho<1/4)
        delta=(1/4)*delta;
    else
        if rho>3/4 && norm(p)== delta
            delta= min(2*delta,1.5);
        else
            delta=delta;
        end
    end
    if (rho>eta)
        xk=xk+p;
    else
        xk=xk;
        
    end
    if norm(g)<tol
        x=xk;
        %disp(x);
        msg="S� se encontr� un m�nimo.";
        bandera=true;
        break;
    end
    k=k+1;
end
if(k==itmax && bandera==false)
    msg="No se encontr� ning�n m�ximo";
    x=xk;
end
end