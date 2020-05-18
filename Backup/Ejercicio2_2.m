%Ejercicio 2.2:
diary('Output_Exercise2_2.txt')
diary on
c1 = 1e-4;
c2 = 0.99;
tol = 1e-5;
eta = 0.1;
DeltaMax = 1.25;
itmax = 1000;
nvalues = [2;8;32;128];
m = length(nvalues);

f = @(x) frosenbrock(x);
    
fprintf("------------Resultados del Problema Rosenbrock------------\n")


%n=2
n = 2;
x0 = ones(n,1);
%Generamos x0
for i=1:n
       if mod(i,2)~=0
           x0(i)=-1.2;
       end
end

%Buscamos el numero de iteraciones totales de cada metodo
[xk, it_SR1] = rcSR1(f, x0,itmax);
[xk, it_BFGS] = lsBFGS( f, x0, itmax);
[xk, it_Mem] = lineLM_BFGS( f, x0,itmax,3);

% Buscamos los valores ||delta(f(xk))|| , f(xk), los errores ||xk - x*|| y 
%el tiempo para las ultimas 5 iteraciones. 

%Para rcSCR1
%vectores auxiliar para las iteraciones
ng_SCR1 = zeros(5,1); 
f_SCR1 = zeros(5,1);
e_SCR1 = zeros(5,1);
t_SCR1 = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);

%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_SR1-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = rcSR1(f,x0,aux(k));
    time = toc;
    ng_SCR1(k) = norm(apGrad(f,xk));
    f_SCR1(k) = f(xk);
    e_SCR1(k) = norm(xk-xaux);
    t_SCR1(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con SCR y n=2 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_SCR1(j), f_SCR1(j), e_SCR1(j), t_SCR1(j))
end

%Para lsBFGS
%vectores auxiliar para las iteraciones
ng_BFGS = zeros(5,1); 
f_BFGS = zeros(5,1);
e_BFGS = zeros(5,1);
t_BFGS = zeros(5,1);
aux = zeros(5,1);
xaux = ones(2,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_BFGS-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lsBFGS(f,x0,aux(k));
    time = toc;
    ng_BFGS(k) = norm(apGrad(f,xk));
    f_BFGS(k) = f(xk);
    e_BFGS(k) = norm(xk-xaux);
    t_BFGS(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsBFGS y n=2 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_BFGS(j), f_BFGS(j), e_BFGS(j), t_BFGS(j))
end

%Para lineLM_BFGS
%vectores auxiliar para las iteraciones
ng_Mem = zeros(5,1); 
f_Mem = zeros(5,1);
e_Mem = zeros(5,1);
t_Mem = zeros(5,1);
aux = zeros(5,1);
xaux = ones(2,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_Mem-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lineLM_BFGS(f,x0,aux(k),3);
    time = toc;
    ng_Mem(k) = norm(apGrad(f,xk));
    f_Mem(k) = f(xk);
    e_Mem(k) = norm(xk-xaux);
    t_Mem(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsLM_BFGS y n=2 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_Mem(j), f_Mem(j), e_Mem(j), t_Mem(j))
end

%Para n=8
n = 8;
x0 = ones(n,1);
%Generamos x0
for i=1:n
       if mod(i,2)~=0
           x0(i)=-1.2;
       end
end
%Buscamos el numero de iteraciones para tomar las ultimas 5 iteraciones
[xk, it_SR1] = rcSR1(f, x0,itmax);
[xk, it_BFGS] = lsBFGS( f, x0, itmax);
[xk, it_Mem] = lineLM_BFGS( f, x0,itmax,3);

% Buscamos los valores ||delta(f(xk))|| , f(xk), los errores ||xk - x*|| y el
% tiempo para las cinco ultimas iteraciones. 

%Para rcSCR1
%vectores auxiliar para las iteraciones
ng_SCR1 = zeros(5,1); 
f_SCR1 = zeros(5,1);
e_SCR1 = zeros(5,1);
t_SCR1 = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_SR1-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = rcSR1(f,x0,aux(k));
    time = toc;
    ng_SCR1(k) = norm(apGrad(f,xk));
    f_SCR1(k) = f(xk);
    e_SCR1(k) = norm(xk-xaux);
    t_SCR1(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con SCR y n=8 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_SCR1(j), f_SCR1(j), e_SCR1(j), t_SCR1(j))
end

%Para lsBFGS
%vectores auxiliar para las iteraciones
ng_BFGS = zeros(5,1); 
f_BFGS = zeros(5,1);
e_BFGS = zeros(5,1);
t_BFGS = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_BFGS-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lsBFGS(f,x0,aux(k));
    time = toc;
    ng_BFGS(k) = norm(apGrad(f,xk));
    f_BFGS(k) = f(xk);
    e_BFGS(k) = norm(xk-xaux);
    t_BFGS(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsBFGS y n=8 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_BFGS(j), f_BFGS(j), e_BFGS(j), t_BFGS(j))
end

%Para lineLM_BFGS
%vectores auxiliar para las iteraciones
ng_Mem = zeros(5,1); 
f_Mem = zeros(5,1);
e_Mem = zeros(5,1);
t_Mem = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_Mem-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lineLM_BFGS(f,x0,aux(k),3);
    time = toc;
    ng_Mem(k) = norm(apGrad(f,xk));
    f_Mem(k) = f(xk);
    e_Mem(k) = norm(xk-xaux);
    t_Mem(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsLM_BFGS y n=8 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_Mem(j), f_Mem(j), e_Mem(j), t_Mem(j))
end

%Para n=32
n = 32;
x0 = ones(n,1);
%Generamos x0
for i=1:n
       if mod(i,2)~=0
           x0(i)=-1.2;
       end
end
%Buscamos el numero de iteraciones para tomar las ultimas 5 iteraciones
[xk, it_SR1] = rcSR1(f, x0,itmax);
[xk, it_BFGS] = lsBFGS( f, x0, itmax);
[xk, it_Mem] = lineLM_BFGS( f, x0,itmax,3 );

% Buscamos los valores ||delta(f(xk))|| , f(xk), los errores ||xk - x*|| y el
% tiempo para las cinco ultimas iteraciones. 

%Para rcSCR1
%vectores auxiliar para las iteraciones
ng_SCR1 = zeros(5,1); 
f_SCR1 = zeros(5,1);
e_SCR1 = zeros(5,1);
t_SCR1 = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_SR1-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = rcSR1(f,x0,aux(k));
    time = toc;
    ng_SCR1(k) = norm(apGrad(f,xk));
    f_SCR1(k) = f(xk);
    e_SCR1(k) = norm(xk-xaux);
    t_SCR1(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con SCR y n=32 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_SCR1(j), f_SCR1(j), e_SCR1(j), t_SCR1(j))
end

%Para lsBFGS
%vectores auxiliar para las iteraciones
ng_BFGS = zeros(5,1); 
f_BFGS = zeros(5,1);
e_BFGS = zeros(5,1);
t_BFGS = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_BFGS-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lsBFGS(f,x0,aux(k));
    time = toc;
    ng_BFGS(k) = norm(apGrad(f,xk));
    f_BFGS(k) = f(xk);
    e_BFGS(k) = norm(xk-xaux);
    t_BFGS(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsBFGS y n=32 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_BFGS(j), f_BFGS(j), e_BFGS(j), t_BFGS(j))
end

%Para lineLM_BFGS
%vectores auxiliar para las iteraciones
ng_Mem = zeros(5,1); 
f_Mem = zeros(5,1);
e_Mem = zeros(5,1);
t_Mem = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_Mem-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lineLM_BFGS(f,x0,aux(k),3);
    time = toc;
    ng_Mem(k) = norm(apGrad(f,xk));
    f_Mem(k) = f(xk);
    e_Mem(k) = norm(xk-xaux);
    t_Mem(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsLM_BFGS y n=32 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_Mem(j), f_Mem(j), e_Mem(j), t_Mem(j))
end

%Para n=128
n = 128;
x0 = ones(n,1);
%Generamos x0
for i=1:n
       if mod(i,2)~=0
           x0(i)=-1.2;
       end
end
%Buscamos el numero de iteraciones para tomar las ultimas 5 iteraciones
[xk, it_SR1] = rcSR1(f, x0,itmax);
[xk, it_BFGS] = lsBFGS( f, x0, itmax);
[xk, it_Mem] = lineLM_BFGS( f, x0,itmax,3);

% Buscamos los valores ||delta(f(xk))|| , f(xk), los errores ||xk - x*|| y el
% tiempo para las cinco ultimas iteraciones. 

%Para rcSCR1
%vectores auxiliar para las iteraciones
ng_SCR1 = zeros(5,1); 
f_SCR1 = zeros(5,1);
e_SCR1 = zeros(5,1);
t_SCR1 = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_SR1-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = rcSR1(f,x0,aux(k));
    time = toc;
    ng_SCR1(k) = norm(apGrad(f,xk));
    f_SCR1(k) = f(xk);
    e_SCR1(k) = norm(xk-xaux);
    t_SCR1(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con SCR y n=128 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_SCR1(j), f_SCR1(j), e_SCR1(j), t_SCR1(j))
end

%Para lsBFGS
%vectores auxiliar para las iteraciones
ng_BFGS = zeros(5,1); 
f_BFGS = zeros(5,1);
e_BFGS = zeros(5,1);
t_BFGS = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_BFGS-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lsBFGS(f,x0,aux(k));
    time = toc;
    ng_BFGS(k) = norm(apGrad(f,xk));
    f_BFGS(k) = f(xk);
    e_BFGS(k) = norm(xk-xaux);
    t_BFGS(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsBFGS y n=128 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_BFGS(j), f_BFGS(j), e_BFGS(j), t_BFGS(j))
end

%Para lineLM_BFGS
%vectores auxiliar para las iteraciones
ng_Mem = zeros(5,1); 
f_Mem = zeros(5,1);
e_Mem = zeros(5,1);
t_Mem = zeros(5,1);
aux = zeros(5,1);
xaux = ones(n,1);
%Creo vector con las ultimas 5 iteraciones
for j = 1:5
    aux(end-j+1) = it_Mem-j+1;
end

k=1;
while k <= length(aux)
    tic;
    [xk,iter] = lineLM_BFGS(f,x0,aux(k),3);
    time = toc;
    ng_Mem(k) = norm(apGrad(f,xk));
    f_Mem(k) = f(xk);
    e_Mem(k) = norm(xk-xaux);
    t_Mem(k) = time;
    k=k+1;
end

%Construimos tabla:
fprintf("------------Resultados con lsLM_BFGS y n=128 ------------\n")
fprintf("\n n \t Iteracion \t ||delta(f(xk))|| \t  f(xk) \t Error \t \t \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \t %d \n", n , aux(j), ng_Mem(j), f_Mem(j), e_Mem(j), t_Mem(j))
end

diary off



    
    
    
    
    

    
    
    











