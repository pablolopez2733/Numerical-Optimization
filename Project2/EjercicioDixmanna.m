close all; clc;
format long;



%% Para n=240
%Create the x0 vector
n=240;
results_iter=zeros(1,5);
results_norm=zeros(1,5);
results_last=zeros(1,5);
results_time=zeros(1,5);
x0 = zeros(n,1);
for k = 1:n
    x0(k) = 2;
end

% CREATE THE DIXMAANA FUNCTION

f=@(x) Dixmanna(x);

vectorm=[1,3,5,17,29];
itmax=5000;
for j=1:5

tic;
[x, iter] = lineLM_BFGS( f, x0, itmax, vectorm(j));
timeElapsed = toc;

results_iter(j)=iter;
results_norm(j)=norm(apGrad(f,x));
results_last(j)=f(x);
results_time(j)=timeElapsed;



end
fprintf('RESULTADOS PARA n=240 : ');
fprintf('\nNumero de iteraciones: ');
fprintf('%d ', results_iter);
fprintf('\nNormas del gradiente: ');
fprintf('%d ', results_norm);
fprintf('\nValores de f(x): ');
fprintf('%d ', results_last);
fprintf('\nTiempos de ejecucion')
fprintf('%d ', results_time);
fprintf('\n');

%% Para n=960:
%Create the x0 vector
n=960;
results_iter2=zeros(1,5);
results_norm2=zeros(1,5);
results_last2=zeros(1,5);
results_time2=zeros(1,5);
x02 = zeros(n,1);
for k = 1:n
    x02(k) = 2;
end

% CREATE THE DIXMAANA FUNCTION

f=@(x) Dixmanna(x);

vectorm2=[1,3,5,17,29];
itmax2=5000;
for j=1:5

tic;
[x2, iter2] = lineLM_BFGS( f, x02, itmax2, vectorm2(j));
timeElapsed = toc;

results_iter2(j)=iter2;
results_norm2(j)=norm(apGrad(f,x2));
results_last2(j)=f(x2);
results_time2(j)=timeElapsed;



end
fprintf('RESULTADOS PARA n=960 : ');
fprintf('\nNumero de iteraciones: ');
fprintf('%d ', results_iter2);
fprintf('\nNormas del gradiente: ');
fprintf('%d ', results_norm2);
fprintf('\nValores de f(x): ');
fprintf('%d ', results_last2);
fprintf('\nTiempos de ejecucion')
fprintf('%d ', results_time2);
fprintf('\n');

%% Tabla n=240
fprintf("------------Resultados con n=240 ------------\n")
fprintf("\n m \t Iteraciones \t ||delta(f(xk))|| \t  f(xk) \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \n", vectorm(j) , results_iter(j), results_norm(j), results_last(j), results_time2(j))
end
%% Tabla n=960
fprintf("------------Resultados con n=960 ------------\n")
fprintf("\n m \t Iteraciones \t ||delta(f(xk))|| \t  f(xk) \t Tiempo \n")
for j = 1:5
fprintf("\n %d \t %d \t %d \t  %d \t %d \n", vectorm(j) , results_iter2(j), results_norm2(j), results_last2(j), results_time(j))
end



