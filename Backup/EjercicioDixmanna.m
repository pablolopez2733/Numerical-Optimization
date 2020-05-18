close all; clc;
format long;



%% Create the x0 vector
n=240;
results_iter=zeros(1,5);
results_norm=zeros(1,5);
results_last=zeros(1,5);
results_time=zeros(1,5);
x0 = zeros(1,n);
for k = 1:n
    x0(k) = 2;
end

%% CREATE THE DIXMAANA FUNCTION

f=@Dixmanna;


%% Ejercicio2.3
vectorm=[1,3,5,17,29];
itmax=20000;
for j=1:5

tic;
[x, iter] = lineLM_BFGS( f, transpose(x0), itmax, vectorm(j));
timeElapsed = toc;

results_iter(j)=iter;
results_norm(j)=norm(apGrad(f,x));
results_last(j)=f(x);
results_time(j)=timeElapsed;

end



