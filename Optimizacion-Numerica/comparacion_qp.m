%Script de comparación
tau = 0.75;
for n = 20:30:1000
    m=floor(n/2);
    p=floor(n/3);
    [Q,A,F,c,b,d]=Generapci(n,m,p,tau);
    t=cputime;
    [x,y,mu,z]=qpinpoint(Q,A,F,c,b,d);
    s=cputime;
    bag1 = [bag; s-t];
end

w=[20:30:1000]';
plot(w,bag1,'dr',w,bag1,'b','Linewidth',3)
title('CPU time en OP')

